from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
import time
import re
import csv
from selenium.webdriver.chrome.options import Options

# 初始化WebDriver
def get_current(x,y,month):
    service = Service(ChromeDriverManager().install())
    chrome_options = Options()
    chrome_options.add_argument("--headless") # 设置无头模式
    
    driver = webdriver.Chrome(service=service,options=chrome_options)
    
    # 打开URL
    url_template = "https://earth.nullschool.net/zh-cn/#2023/0"+str(month)+"/01/0000Z/ocean/surface/currents/anim=off/orthographic=12.05,38.38,5382/loc="
    url = url_template + str(x) + "," + str(y)
    
    driver.get(url)

    # 等待JavaScript加载完成，可能需要调整等待时间
    driver.implicitly_wait(20) # 10秒

    time.sleep(3)
    # 使用By.XPATH来定位元素
    element = driver.find_element(By.XPATH, '//*[@id="spotlight-panel"]/div[2]/div')
    print(element.get_attribute('innerHTML'))
    content = element.get_attribute('aria-label')
    matches = re.findall(r"(\d+°) @ ([\d.]+)", content)
    driver.quit()
    if matches:
        degree, value = matches[0] # 提取第一组匹配结果
        degree = int(degree.replace("°", ""))
        value = float(value)
    else:
        degree, value = None, None
        

# 将值转换为浮点数
    
    
    
    return degree, value
# 关闭浏览器
# file_path = 'current_April.csv'
# x_list = [i/100 for i in range(int(18.57*100), int(20*100) + 1, 10)]#15.28
# y_list = [i/100 for i in range(int(35.82*100), int(40.52*100) + 1, 10)]
# for x in x_list:
#     for y in y_list:
#         with open(file_path, mode='a', newline='') as file:
#             direction, velocity = get_current(x,y,4)
#             data_to_append = [x, y, direction, velocity]
#             writer = csv.writer(file)
#             writer.writerow(data_to_append)
            
file_path = 'current_July.csv'
x_list = [i/100 for i in range(int(17.17*100), int(20*100) + 1, 10)]#15.28
y_list = [i/100 for i in range(int(35.82*100), int(40.52*100) + 1, 10)]
for x in x_list:
    for y in y_list:
        with open(file_path, mode='a', newline='') as file:
            direction, velocity = get_current(x,y,7)
            data_to_append = [x, y, direction, velocity]
            writer = csv.writer(file)
            writer.writerow(data_to_append)
file_path = 'current_Oct.csv'
for x in x_list:
    for y in y_list:
        with open(file_path, mode='a', newline='') as file:
            direction, velocity = get_current(x,y,10)
            data_to_append = [x, y, direction, velocity]
            writer = csv.writer(file)
            writer.writerow(data_to_append)