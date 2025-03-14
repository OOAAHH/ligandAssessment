# -*- coding: utf-8 -*-
# input ligand_mass_center.json
# output grid unit distance
# Author: HaoSun
# Date: 2025-03-11 15:24
import os
import json
import numpy as np
import pandas as pd
import math

def convert_to_grid_unit(coordinate, spacing=2.0):
    """
    使用新的逻辑将单个坐标转换为2Å网格单元坐标。
    对于 spacing=2，网格单元 n 对应真实坐标区间 [2*n - 1, 2*n + 1)
    因此采用公式： n = floor((c + spacing/2) / spacing)
    """
    return [math.floor((c + spacing/2) / spacing) for c in coordinate]

def calculate_distance(grid_unit_ref, grid_unit_predict):
    """
    计算参考结构和预测结构的 grid unit 坐标的距离
    计算公式：|x_ref - x_predict| + |y_ref - y_predict| + |z_ref - z_predict|
    """
    distance = sum(abs(np.array(grid_unit_ref) - np.array(grid_unit_predict)))
    return distance

def read_config(config_file):
    """
    读取配置文件，返回参考结构名称和ligand选择条件。
    
    配置文件要求：
      第一行：参考结构名称（不带扩展名），用于区分参考与predict结构
      第二行：ligand选择条件（例如 "chain A" 或 "resn LIG"）
    """
    with open(config_file, 'r') as f:
        lines = f.readlines()
    if len(lines) < 2:
        raise ValueError("配置文件至少需要两行：第一行参考结构名称，第二行ligand选择条件")
    reference_name = lines[0].strip()  # 参考结构名称
    select_ligand_chain = lines[1].strip()  # ligand选择条件
    return reference_name, select_ligand_chain

def process_ligand_json(input_folder, config_file):
    """
    读取输入的 JSON 文件，转换坐标为 grid unit，并计算每个预测结构与参考结构之间的距离。
    返回一个包含转换后的坐标和距离的表格数据。
    """
    # 检查cache文件夹是否存在
    cache_folder = os.path.join(input_folder, "cache")
    if not os.path.exists(cache_folder):
        print("Cache folder does not exist. Please run massC calculation first.")
        return [], None, cache_folder  # 如果cache文件夹不存在，返回None
    
    # 检查缓存中的ligand质心文件
    input_json = os.path.join(cache_folder, "ligand_mass_center.json")
    
    if not os.path.exists(input_json):
        print(f"ligand_mass_center.json not found in {cache_folder}. Please run massC calculation first.")
        return [], None, cache_folder  # 如果没有文件，则返回None

    # 解析 JSON 文件
    with open(input_json, 'r') as f:
        data = json.load(f)
    
    # 从配置文件中读取参考结构名称
    reference_name, _ = read_config(config_file)

    # 获取参考结构的坐标
    if reference_name not in data:
        print(f"Reference structure {reference_name} not found in the JSON data.")
        return [], None, cache_folder  # 如果参考结构不在数据中，返回空列表和None
    
    reference_coords = data[reference_name]  # 参考结构的坐标

    # 计算参考结构的 grid unit 坐标
    reference_grid_unit = convert_to_grid_unit(reference_coords)

    # 准备表格数据
    table_data = []

    # 计算每个预测结构的 grid unit 坐标，并与参考结构进行比较
    for name, coords in data.items():
        if name == reference_name:  # 排除参考结构
            continue

        grid_unit = convert_to_grid_unit(coords)
        distance = calculate_distance(reference_grid_unit, grid_unit)
        
        # 为每个结构名称记录坐标和距离
        table_data.append([name, grid_unit, distance])

    return table_data, reference_grid_unit, cache_folder

def save_to_csv(table_data, reference_grid_unit, output_file):
    """
    将表格数据保存到 CSV 文件。
    """
    # 添加表头
    columns = ["Structure Name", "Grid Unit", "Distance"]

    # 将数据转换为 DataFrame 并保存为 CSV
    df = pd.DataFrame(table_data, columns=columns)
    
    # 添加参考结构的 Grid Unit 坐标作为表头的一部分
    reference_row = ["Reference", reference_grid_unit, "NA"]
    df.loc[-1] = reference_row  # 添加参考结构行
    df.index = df.index + 1  # 调整索引
    df.sort_index(inplace=True)

    # 保存 CSV 文件
    df.to_csv(output_file, index=False)

def find_ligan_massc_grid_units(input_folder, config_file):
    """
    主函数，读取 JSON 数据，计算每个预测结构的 Grid Unit 坐标与参考结构的距离，
    并将结果保存到 CSV 文件。
    """
    # 处理输入 JSON 文件
    table_data, reference_grid_unit, cache_folder = process_ligand_json(input_folder, config_file)

    if not table_data:
        return  # 如果没有处理数据，直接返回

    # 构造输出文件路径
    output_file = os.path.join(cache_folder, "ligand_massC_grid_unit_comparison_with_distance.csv")
    
    # 保存结果到 CSV 文件
    save_to_csv(table_data, reference_grid_unit, output_file)
    print(f"Results have been saved to {output_file}")
