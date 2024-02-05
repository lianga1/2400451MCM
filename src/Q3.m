% 已知一个轨迹
X = readtable('path.csv');
X.Properties.VariableNames = {'x','y','z','dx','dy','dz'};

% 取点,得到每个点的先验概率
% 根据表格中的数据，均匀取得6个点，存储在P点中
% Assuming 'yourTable' is your n x 6 table
nRows = size(X, 1); % Get the number of rows in the table

% Generate 6 unique row indices for sampling
sampleIndices = sort(randperm(nRows, 6));

% Sample the rows
% Gaussian parameters
mu_x = X{end, 'x'}; sigma_x = 2000;
mu_y = X{end, 'y'}; sigma_y = 2000;
mu_z = X{end, 'z'}; sigma_z = 200;

% Function to calculate Gaussian PDF
gaussianPDF = @(x, mu, sigma) (1/(sigma*sqrt(2*pi))) * exp(-((x-mu).^2)/(2*sigma^2));

% Sample points
samples = X{sampleIndices, :};

% Calculate probabilities
% 假设探测范围为搜索点为原点方球1000米。
probabilities = zeros(1,size(samples, 1));
for i = 1:size(samples, 1)
    % Integrate the Gaussian PDF over the specified range for each dimension
    prob_x = integral(@(x) gaussianPDF(x, mu_x, sigma_x), samples(i,1)-500, samples(i,1)+500);
    prob_y = integral(@(y) gaussianPDF(y, mu_y, sigma_y), samples(i,2)-500, samples(i,2)+500);
    prob_z = integral(@(z) gaussianPDF(z, mu_z, sigma_z), -500, 0);
    
    % Multiply the probabilities to get the combined probability for each point
    probabilities(i) = prob_x * prob_y * prob_z;
end

p_search = 0.95;
k = 0.0011;
exp_prob = @(p,r) p*exp(-k*r)*k; 
% 遍历每个点，进行随机数搜寻，把每个点的概率变化记录下来
all_probs = [];
P = zeros(1,size(probabilities,2));
accumulate_prob = [];
times = [];
for i =1:1000
    accumulate_prob_temp = zeros(1,size(probabilities,2));
    for j=1:size(probabilities,2)
        %生成0，1之间的均匀随机数n
        n = rand;
        P(j) = probabilities(j)*integral(@(r) exp_prob(p_search,r),0,1000)^3;
        
        if n < P
            disp(['searched']);
            times = [times,i];
            disp(i);
            
            
        else
            prob_temp = probabilities(j);
            for k = 1:size(probabilities,2)
                if k == j
                    probabilities(k) = prob_temp/6.01;
                else
                    probabilities(k) = probabilities(k) + prob_temp/6.01; 
                end
            end
        end
        for a = 1:size(accumulate_prob_temp,2)
            if size(all_probs,2) >0
            accumulate_prob_temp(a) = sum(all_probs(1:a:end,a))+P(a);
            else
                accumulate_prob_temp(a) = P(a);
            end
            
        end

    end
    all_probs = [all_probs;P];
    accumulate_prob = [accumulate_prob;accumulate_prob_temp];
end
%画图
accumulate_prob = accumulate_prob/9;
%把所有的all_probs，画成图，每一列为一个点的概率变化，每一行为一个时间点
% Assuming 'data' is your 100x10 matrix

% Plot each column as a separate line
figure;
hold on;
colors = lines(10); % Generate 10 distinct colors
for i = 1:size(all_probs, 2)
    plot(1:1000, all_probs(:, i), 'Color', colors(i, :), 'LineWidth', 2);
end
for i =1:5
xline(times(i), '-'); % 为每个数据点绘制一条竖直虚线
end
hold off;

% Enhancing the plot
xlabel('Time Step');
ylabel('Value');
title('Evolution of 6 Points Over Time');
legend('Point 1', 'Point 2', 'Point 3', 'Point 4', 'Point 5', 'Point 6');% ...
    %'Point 7', 'Point 8', 'Point 9', 'Point 10');

figure;
hold on;
colors = lines(10); % Generate 10 distinct colors
for i = 1:size(accumulate_prob, 2)
    plot(1:1000, accumulate_prob(:, i), 'Color', colors(i, :), 'LineWidth', 2);
end
for i = 1:5
xline(times(i), '-'); % 为每个数据点绘制一条竖直虚线
end
hold off;

% Enhancing the plot
xlabel('Time Step');
ylabel('Value');
title('accumulate of 6 Points Over Time');
legend('Point 1', 'Point 2', 'Point 3', 'Point 4', 'Point 5', 'Point 6');% ...
    %'Point 7', 'Point 8', 'Point 9', 'Point 10');
