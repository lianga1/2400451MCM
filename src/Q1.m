rho = 1026 ;%kg/m3
mu = 1e-3; %粘性系数
m = 20000 ;%kg
V = 40 ;%m^3
m_w = 20000;%kg
global velocity;
global depths;
global direction;
global zone;
direction = [0 0 pi];
%定义读取洋流速度的模型

velocity = get_current_velocity('current.csv');
depths = get_depth('depth.csv');
KalmanFilter();
function KalmanFilter()
    % 初始化状态和协方差
    coord = [18.07,37.92,-200];
    xy = [1,2];
    [xy(1),xy(2)]= coord2x(coord(1),coord(2));
    x_hat = [xy(1) xy(2) -200 3 5 -0.2];
    %x_hat = [18.07; 37.92; 200; 0; 0; 0]; % 初始状态估计，包含经纬度、深度和初始速度
    P = [1.1,0,0,0,0,0;
         0.1,1.1,0,0,0,0;
         0,0,1.1,0,0,0;
         0,0,0,1.1,0.1,0;
         0,0,0,0.1,1.1,0;
         0,0,0,0,0,1;]; % 初始误差协方差矩阵
    Q = diag([0.1 0.1 0.1 0.1 0.1 0.1]); % 过程噪声协方差矩阵
    R = diag([10 10 10]); % 观测噪声协方差矩阵，假设只观测到位置和深度
    
    % 控制输入（如果有的话）
    u = []; % 此处留空，因为您的模型似乎没有直接的控制输入
    H = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 0 0];
    all_data=[];
    % 获取观测数据（假设已经有了某种形式的观测数据）
    X2_all = [];
    Z2_all = [];
    all_coord = [];
    for k = 1:100
        % 预测步骤
        [~,X] = solve(@ode, x_hat); % 使用ODE求解预测下一状态
        data = [];
        for i = 1:size(X,1)
            data_temp=zeros(1,6);
            [data_temp(1),data_temp(2)]=x2coord(X(i,1),X(i,2));
            data_temp(3)=X(i,3);
            data = [data;data_temp];
        end
        all_coord = [all_coord;data];
        x_hat_1 = x_hat;
        [~,X2 ]=solve(@ode,x_hat_1);
        for i = 1:size(X2,2)
        X2(i,:) =getObservations(X2(i,:));
        
        end
        X2_all = [X2_all;X2];
        x_hat_minus = X(end,:); % 获取ODE求解的最终状态作为预测状态
        
        Z = getObservations(x_hat(1:6)); % 假设这个函数返回观测数据
        % 预测误差协方差（这里简化处理，实际应用中可能需要更复杂的模型来估计P）
        P_minus = P + Q;
        
        % 更新步骤
        K = P_minus * H' / (H * P_minus * H' + R); % 计算卡尔曼增益，需要定义H
        z = [Z(1);Z(2);Z(3)]; % 当前时刻的观测数据
        z2=[Z(1),Z(2),Z(3)];
        Z2_all = [Z2_all;z2];
        x_hat = x_hat_minus' + K * (z - H * x_hat_minus'); % 更新状态估计
        P = (eye(6) - K * H) * P_minus; % 更新误差协方差
        disp(['Time Step ', num2str(k), ':']);
        disp(['Predicted State: ', num2str(x_hat')]);
        all_data = [all_data;X];
        if k == 34
            disp('pause');
        end
    end
    all_data = all_data -0.2;
    figure;
    hold on;
    scatter(all_data(1:2:end,1), all_data(1:2:end,2), 'r', 'DisplayName', 'Kalman'); % Red for Table 1
    scatter(X2_all(2:2:end,1), X2_all(2:2:end,2),'Color',[1, 0.7529, 0], 'DisplayName', 'Predict'); % Green for Table 2
    scatter(Z2_all(:,1), Z2_all(:,2), 'Color',[70/255, 114/255, 196/255], 'DisplayName', 'Observe'); % Blue for Table 3
    hold off;

    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    title('Scatter Plots of Kalman Filter');
    legend;
    figure;
    hold on;
    plot(all_data(1:2:end,1), all_data(1:2:end,2), 'r-', 'LineWidth', 2); % Plot table1 in red
    plot(X2_all(2:2:end,1), X2_all(2:2:end,2) ,'Color',[1, 0.7529, 0], 'LineWidth', 2); % Plot table2 in green
    plot(Z2_all(:,1), Z2_all(:,2), 'Color',[70/255, 114/255, 196/255], 'LineWidth', 2); % Plot table3 in blue
    hold off;
    
    % Enhancements
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    title('Route of three');
    legend('Kalman', 'Predict', 'Observe');
    writematrix(all_coord,'path.csv');

end

function Z = getObservations(x_hat)
    % 获取或生成观测数据
    sigma = sqrt(0.1);
    
    % 生成正态分布的随机扰动
    noise = sigma * randn(1, size(x_hat,6)); % 生成三个随机数，均值为0，标准差为sigma
    
    % 生成新的向量z
    Z = x_hat + noise;
    
end


function [x,y]=coord2x(lat,lon)
    global zone;
    zone = utmzone(lat, lon);
    
    % 配置UTM投影参数
    utmstruct = defaultm('utm');
    utmstruct.zone = zone;
    utmstruct.geoid = wgs84Ellipsoid(); % 使用WGS84椭球体
    utmstruct = defaultm(utmstruct);
    
    % 将经纬度转换为UTM坐标
    [x,y] = mfwdtran(utmstruct, lat, lon);
end

function [coord_x ,coord_y]=x2coord(x,y)
    global zone;
    
    utmstruct = defaultm('utm');
    utmstruct.zone = zone;
    utmstruct.geoid = wgs84Ellipsoid(); % 使用WGS84椭球体
    utmstruct = defaultm(utmstruct);
    [coord_x,coord_y] = minvtran(utmstruct, x, y);
end


function velocity=get_current_velocity(filename)
    % 假设您的CSV文件名为"data.csv"
    data = readtable(filename);
    data.Properties.VariableNames = {'x', 'y', 'degree', 'speed'};
    % 使用readtable函数读取文件

    % 检查并替换所有NaN值为指定的占位符，例如用-999表示空值
    placeholder = -999;
    velocity = fillmissing(data, 'constant', placeholder);
    
    % 显示前几行数据以验证
    

end

function [t,X] = solve(ode, init)
% SOLVE 解三维矢量的一元二阶常微分方程
%   ode: 一个函数句柄，表示ODE右侧的f(x, dx/dt, t)，这里x是三维矢量
%   init: 初始条件向量，[x0 dx0 y0 dy0 z0 dz0]形式

% 将二阶ODE转换为一阶ODE系统的函数
function dydt = odeSystem(t, y)
    dydt = zeros(6,1); % 初始化输出为6x1向量
    f = ode(y(1:3), y(4:6), t); % 调用用户提供的函数计算加速度
    dydt(1:3) = y(4:6); % 更新位置的导数
    dydt(4:6) = f; % 更新速度的导数
end

% 初始条件
y0 = init; % 这里假设init已经是[ x0 y0 z0 dx0 dy0 dz0 ]格式

% 时间跨度
tspan = [0 10]; % 根据需要调整

% 使用ode45求解
[t,X] = ode45(@odeSystem, tspan, y0);

% 结果
%xt = [t y]; % 将时间和解（x, y, z的值）组合在一起返回
end

function a = ode(x, dxdt, t)
    % 示例ODE函数，计算三维空间中一个物体的加速度
    % x: 位置向量，格式为[x; y; z]
    % dxdt: 速度向量，格式为[dx/dt; dy/dt; dz/dt]
    % t: 时间，本示例中未使用
    % a: 加速度向量，格式为[ax; ay; az]

    % 假设的恒定力F和物体质量m
    g= [0; 0; -9.8]; % 仅在z方向上有重力作用，模拟简单下落
    m = 20000 ;%kg
    rho = 1026 ;%kg/m3
    mu = 1e-3; %粘性系数
    V = 40 ;%m^3
    m_w = 20000;%kg
    C=-0.47;
    l=15;b=3;c=2;
    s=getsquare(l,b,c);
    [x1,y1] = x2coord(x(1),x(2));
    coord = [x1,y1,x(3)];
    v_c = getcurrent(coord);
    %if dxdt < 1
    %    F_p=[60000;0;0];
    %elseif dxdt <2
    %    F_p = [100000;0;0];
    %else
    F_p=[367749;0;0];
    %end

    F_m = F_p ./dxdt;
    if dxdt(3)==0
        F_m(3) = 0;
    end
    if v_c == -999
        return
    end
    F = ((m+m_w).*g - rho*V.*g)*0 +0.5*C*rho*dxdt.^2.*s'.*sign(dxdt) + 0.5*C*rho*v_c'.^2.*s'.*sign(v_c')+F_m;
    % 根据牛顿第二定律计算加速度
    if dxdt(1)<3
        F(1) = m+m_w;
    end
    a = F / (m+m_w);
    a(3)=0;

    % 本示例中，加速度a与位置x和速度dxdt无关，仅为演示
    % 在更复杂的情况下，a可能依赖于x, dxdt和t
end

function s = getsquare(a,b,c)%宽b，长a，高c
    %x方向上的速度
    s_x = 0.5*pi*b*c;
    s_y = 0.5*pi*a*c;
    s_z = 0.5*pi*a*b;
    s = [s_x s_y s_z];
end

function v_c = getcurrent(x)
% 根据深度层流模型和现有数据，进行建模
    %模型为$\frac{dP}{dx}=\mu \frac{d^2v}{dy^2}$
    global velocity;
    global direction;
    global depths;

    k = 1;
    %fprintf("x %f, y %f\n",x(1),x(2))
    condition = (abs(velocity.x-x(1))<0.09999) & (abs(velocity.y-x(2))<0.09999);
    deg = velocity.degree(condition);
    deg = deg(1);
    spd = velocity.speed(condition);
    spd = spd(1);
    condition_d = abs(depths.x- x(1))<0.09999 & abs(depths.y-x(2))<0.09999;
    depth = depths.depth(condition_d);
    depth = depth(1);
    v_c_surf = -[spd*cos(direction(3)-rad2deg(deg)) spd*sin(direction(3)-rad2deg(deg)) 0];
    if v_c_surf == -999
        fprintf('collapse at point%f,%f\n',x(1),x(2));
        v_c = -999;
        return
    end
    if depth == 0
        fprintf('collapse at point%f,%f\n',x(1),x(2));
        v_c = -999;
        return 
    end
    v_c = v_c_surf * exp(-k * (x(3)/depth));
end

function depths = get_depth(filename)
    data = readtable(filename);
    data.Properties.VariableNames =  {'y', 'x', 'depth'};
    placeholder = 0;
    depths = fillmissing(data, 'constant', placeholder);
    depths = depths;
end

