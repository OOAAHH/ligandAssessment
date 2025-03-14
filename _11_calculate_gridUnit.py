# -*- coding: utf-8 -*-
# output grid unit
# Author: HaoSun
# Date: 2025-03-08 10:22
import math

def convert_to_grid_unit(coordinate):
    """
    将原始坐标转换为 2Å 网格单元坐标。

    原则：左闭右开区间 [0, 2)，坐标将转换为每2Å一个单位。
    对于坐标 (x, y, z)，将其每个分量转换为新的 grid unit 坐标：
    grid_x = ceil(x / 2)  # 对于正数
    grid_y = ceil(y / 2)  # 对于正数
    grid_z = ceil(z / 2)  # 对于正数
    对于负数需要使用 math.floor 来取整。

    参数：
    coordinate : dict
        包含原始坐标的数据，形如 { "id": [x, y, z] }

    返回：
    dict : 转换后的 grid unit 坐标，形如 { "id_gridUnit": [grid_x, grid_y, grid_z] }
    # 示例用法
    coordinate = {
        "6GZR.1-10_02_super": [
            16.99809881148614,
            -0.33309801100142905,
            4.7646988864380315
        ]
    }

    converted_coordinates = convert_to_grid_unit(coordinate)
    print(converted_coordinates)
    """
    converted_coordinates = {}
    
    for key, value in coordinate.items():
        # 对每个坐标分量进行转换
        grid_x = math.floor(value[0] / 2) if value[0] < 0 else math.ceil(value[0] / 2)
        grid_y = math.floor(value[1] / 2) if value[1] < 0 else math.ceil(value[1] / 2)
        grid_z = math.floor(value[2] / 2) if value[2] < 0 else math.ceil(value[2] / 2)

        # 存储转换后的坐标
        new_key = f"{key}_gridUnit"
        converted_coordinates[new_key] = [grid_x, grid_y, grid_z]
    
    return converted_coordinates


