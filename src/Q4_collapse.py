from matplotlib import pyplot as plt
import numpy as np

v_x = 5  # x轴速度, m/s
s_y = 15 * np.pi  # y轴截面积
initial_distance_x = 200  # 初始x轴距离, m
avoidance_distance_x = 30  # 避让距离, m

# 计算避让开始的时间
time_to_avoidance = (initial_distance_x - avoidance_distance_x) / (2 * v_x)

# 椭球体1和椭球体2在避让后2s的时间点
time_post_avoidance = time_to_avoidance + 2

time_continue = time_post_avoidance + 6

# 时间范围
t = np.linspace(0, time_continue, 1000)

# x轴上的位置
x1 = 5 * t
x2 = 200 - 5 * t

# y轴上的位置，在避让开始后
y1 = np.piecewise(t, [t < time_to_avoidance, t >= time_to_avoidance],
                  [lambda t: 0, lambda t: -0.75 * (t - time_to_avoidance)**2])
y2 = np.piecewise(t, [t < time_to_avoidance, t >= time_to_avoidance],
                  [lambda t: 0, lambda t: 0.75 * (t - time_to_avoidance)**2])

# 为简化，2秒后y轴的加速度近似为一个定值
# 假定加速度为正负常数值
a_y1 = 1.5  # 椭球体1的y轴加速度, m/s^2，假设值
a_y2 = -1.5  # 椭球体2的y轴加速度, m/s^2，假设值
# 在2秒后y轴的位置更新
# y(t) = y0 + v0t + 0.5at^2
# y1[t >= time_post_avoidance] = y1[t==time_post_avoidance] + 3 * (t[t >= time_post_avoidance] - time_post_avoidance) + 0.5 * a_y1 * (t[t >= time_post_avoidance] - time_post_avoidance)**2
# y2[t >= time_post_avoidance] = y2[t==time_post_avoidance] + 3 * (t[t >= time_post_avoidance] - time_post_avoidance) + 0.5 * a_y2 * (t[t >= time_post_avoidance] - time_post_avoidance)**2
# y1 = np.piecewise(t,[t < time_post_avoidance, t >= time_post_avoidance],
                #   [lambda t: 0, lambda t: -0.75 * (time_to_avoidance - time_post_avoidance)**2 + -3 * (t - time_post_avoidance) + 0.5 * a_y1 * (t - time_post_avoidance)**2])
# y2 = np.piecewise(t,[t < time_post_avoidance, t >= time_post_avoidance],

                    # [lambda t: 0, lambda t: 0.75 * (time_to_avoidance - time_post_avoidance)**2 + 3 * (t - time_post_avoidance) + 0.5 * a_y2 * (t - time_post_avoidance)**2])
                    # 计算在 time_post_avoidance 之后的 y1 和 y2 的值
y1_start = -0.75 * 2**2  # 计算 time_post_avoidance 时刻的 y1 初始值
y2_start = 0.75 * 2**2  # 计算 time_post_avoidance 时刻的 y2 初始值

# 更新 y1 和 y2 的值，确保形状一致
delta_t = t[t >= time_post_avoidance] - time_post_avoidance  # 从避让开始后到当前的时间差
y1_update = y1_start + -3 * delta_t + 0.5 * a_y1 * delta_t**2
y2_update = y2_start + 3 * delta_t + 0.5 * a_y2 * delta_t**2

y1[t >= time_post_avoidance] = y1_update
y2[t >= time_post_avoidance] = y2_update
y1[t >= time_post_avoidance +2] = -6
y2[t >= time_post_avoidance +2] = 6

# y1[t<=time_continue]=y1_update[-1]
# y2[t<=time_continue]=y2_update[-1]

# 绘制图像
plt.figure(figsize=(12, 6))
plt.plot(x1[t <= time_continue], y1[t <= time_continue], label='Submersible 1')
plt.plot(x2[t <= time_continue], y2[t <= time_continue], label='Submersible 2')
plt.xlabel('X Position (m)')
plt.ylabel('Y Position (m)')
plt.title('Motion Trajectories of Two Collapsing Submersibles')
plt.legend()
plt.grid(True)
plt.show()