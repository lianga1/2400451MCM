# 解决一个旅行商问题
# 问题描述：一个旅行商要拜访n个城市，他要找到一条路径，使得他经过每个城市一次，且最后回到出发的城市。
# 旅行商问题是一个NP难问题，没有多项式时间的解法。但是可以使用动态规划来解决规模较小的问题。

# 问题：给定n个城市之间的距离，求旅行商的最短路径。
# 输入：n个城市之间的距离矩阵，n
# 输出：旅行商的最短路径

def nearest_neighbor_tsp(distance_matrix):
    n = len(distance_matrix)
    visited = [False] * n
    path = [0]  # 假设从城市0开始
    visited[0] = True
    total_distance = 0

    for _ in range(1, n):
        last = path[-1]
        next_city = None
        min_distance = float('inf')
        for i in range(n):
            if not visited[i] and distance_matrix[last][i] < min_distance:
                min_distance = distance_matrix[last][i]
                next_city = i
        path.append(next_city)
        visited[next_city] = True
        total_distance += min_distance

    # 添加从最后一个城市返回到起始城市的距离
    total_distance += distance_matrix[path[-1]][path[0]]
    return total_distance, path

# 使用相同的城市之间的距离矩阵进行测试

# 点1: (51, 92)
# 点2: (14, 71)
# 点3: (60, 20)
# 点4: (82, 86)
# 点5: (74, 74)
distance_matrix = [[  0.        ,  61.40032573,  55.07267925,  29.06888371,  61.40032573,  24.69817807],
 [ 61.40032573,   0.        , 116.21101497,  81.0246876 ,  86.83317338,  84.89994111],
 [ 55.07267925, 116.21101497,   0.        ,  50.03998401,  87.80091116,  35.51056181],
 [ 29.06888371,  81.0246876 ,  50.03998401,   0.        ,  40.31128874,  17.        ],
 [ 61.40032573,  86.83317338,  87.80091116,  40.31128874,   0.        ,  57.30619513],
 [ 24.69817807,  84.89994111,  35.51056181,  17.        ,  57.30619513,   0.        ]]



# 计算最近邻居算法的结果
total_distance, path = nearest_neighbor_tsp(distance_matrix)
print("Total Distance:", total_distance)
print("Path:", path)
