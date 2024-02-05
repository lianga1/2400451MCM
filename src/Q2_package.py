import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler  

file_path = 'quienment data.xlsx'  
df = pd.read_excel(file_path, engine='openpyxl')  

scaler = MinMaxScaler()  
  
# 对每一列进行归一化  
for column in df.columns:  
    df[column] = scaler.fit_transform(df[[column]])  

df["costs"] = (df["purchase"] + df["maintenance"] + df["prepare"] + df["using"]) * 0.8 + df["weight"] * 0.2
df["efficiencies"] = df["scope"] * 0.4 + df["accuary"] * 0.4 + df["failure"] * 0.2  

print(df)



def knapsack(items_matrix, max_cost_factor):
    n = len(items_matrix)  # 物品数量
    W = max_cost_factor  # 背包最大容量
    dp = [[0 for _ in range(W + 1)] for _ in range(n + 1)]  # 初始化动态规划表

    # 动态规划填表
    for i in range(1, n + 1):
        cost = items_matrix[i-1][0]  # 当前物品的成本
        efficiency = items_matrix[i-1][1]  # 当前物品的效率
        for w in range(1, W + 1):
            if cost <= w:
                dp[i][w] = max(dp[i-1][w], dp[i-1][w-cost] + efficiency)
            else:
                dp[i][w] = dp[i-1][w]

    # 回溯找出选择的物品
    max_efficiency = dp[n][W]
    w = W
    selected_items = [False] * n
    for i in range(n, 0, -1):
        if dp[i][w] != dp[i-1][w]:
            selected_items[i-1] = True
            w -= items_matrix[i-1][0]

    return selected_items, max_efficiency

def plan_load_dp(items, max_cost):
    # 将内存和时间转换为整数
    total_cost = sum(items['costs'])
    items['costs'] = items['costs'] /total_cost
    min_cost = min(items['costs'])
    cost_scale = 100 / min_cost

    scaled_cost_effi = [(int(items['costs'][i] * cost_scale), items['efficiencies'][i]) for i in range(len(items['costs']))]
    max_cost_scaled = int(max_cost * cost_scale)

    return knapsack(scaled_cost_effi, max_cost_scaled)
    # 初始化动态规划表
    # n = len(items)
    # dp = [[float(0)] * (max_cost_scaled + 1) for _ in range(n + 1)]
    # dp[0][0] = 0

    # # 动态规划填表
    # for i in range(1, n + 1):
    #     cost, efficiency = scaled_cost_effi[i - 1]
    #     for j in range(max_cost_scaled + 1):
    #         dp[i][j] = max(dp[i][j], dp[i-1][j])
    #         if j >= cost:
    #             dp[i][j] = max(dp[i][j], dp[i-1][j-cost] + efficiency)

    # # 找出最小计算时间增加量对应的内存节约
    # max_efficiency = max(dp[n])
    # for j in range(max_cost_scaled + 1):
    #     if dp[n][j] == max_efficiency:
    #         final_cost = j
    #         break

    # # 通过回溯找出具体的卸载决策
    # load_decision = [False] * n
    # for i in range(n, 0, -1):
    #     if final_cost >= scaled_cost_effi[i - 1][0] and dp[i][final_cost] == dp[i-1][final_cost - scaled_cost_effi[i - 1][0]] + scaled_cost_effi[i - 1][1]:
    #         load_decision[i - 1] = True


    #return load_decision, max_efficiency

# print(plan_load_dp(df, 0.3))
result =[]

for i in range(1, 11):
    di = df.copy(deep=True)
    plan, efficiency = plan_load_dp(di, i/10)
    print(plan, efficiency)
    result.append((i/10, efficiency,plan))

# write the result to a excel file
dr = pd.DataFrame(result, columns=['max_cost', 'efficiency', 'plan'])
dr.to_excel('result.xlsx', index=False, header=True)
    