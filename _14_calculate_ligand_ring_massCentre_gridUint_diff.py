# -*- coding: utf-8 -*-
# output ring grid unit distance
# 读取 ligand_ring_mass_center.json 中的环重原子质心数据，
# 输出各结构中每个环的 grid unit 坐标及与参考结构的距离。
# Author: HaoSun
# Date: 2025-03-12  (更新于新版逻辑)

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

def read_config(config_file):
    """
    读取配置文件，返回参考结构名称。
    配置文件要求：
      第一行：参考结构名称（不带扩展名）
    """
    with open(config_file, 'r') as f:
        lines = f.readlines()
    if len(lines) < 1:
        raise ValueError("配置文件至少需要一行：参考结构名称")
    reference_name = lines[0].strip()
    return reference_name

def calculate_distance(grid_unit_ref, grid_unit_predict):
    """
    计算参考结构和预测结构的 grid unit 坐标的距离
    计算公式：|x_ref - x_predict| + |y_ref - y_predict| + |z_ref - z_predict|
    """
    return sum(abs(np.array(grid_unit_ref) - np.array(grid_unit_predict)))

def find_ligand_ring_grid_units(input_folder, config_file, spacing=2.0):
    """
    读取 ligand_ring_mass_center.json 中存储的环的重原子质心数据，
    对每个结构的每个环计算 grid unit 坐标，并与参考结构进行比较计算距离。
    
    以每个环的 "atoms" 列表拼接成的字符串作为唯一标识符 (ring_id)。
    
    最后构造一个表格：第一列为 ring_id，第二列为参考结构的 grid unit，
    随后依次为每个预测结构的 grid unit 及其与参考结构的距离。
    表格保存为 CSV 文件，文件名为 "ligand_ring_grid_unit_comparison_with_distance.csv"，
    并返回各结构的环 grid unit 坐标字典。
    """
    cache_folder = os.path.join(input_folder, "cache")
    input_json = os.path.join(cache_folder, "ligand_ring_mass_center.json")
    
    if not os.path.exists(input_json):
        print(f"{input_json} 不存在，请先运行 ligand_ring_mass_center 计算。")
        return {}
    
    # 读取 JSON 数据
    with open(input_json, 'r') as f:
        data = json.load(f)
    
    reference_name = read_config(config_file)
    
    # 构建结果字典：结构名称 -> { ring_id : grid unit }
    results = {}
    # 收集所有 ring_id 用于构造表格行
    all_ring_ids = set()
    
    for struct_name, rings in data.items():
        ring_dict = {}
        for ring in rings:
            # 以 atoms 列表拼接成字符串作为 ring_id
            ring_id = "_".join(ring["atoms"])
            all_ring_ids.add(ring_id)
            center = ring["center_of_mass"]
            grid_unit = convert_to_grid_unit(center, spacing)
            ring_dict[ring_id] = grid_unit
        results[struct_name] = ring_dict
    
    # 检查参考结构是否存在
    if reference_name not in results:
        print(f"参考结构 {reference_name} 不存在于 JSON 数据中。")
        return results
    
    reference_rings = results[reference_name]
    
    # 构造表格数据：每行对应一个 ring_id
    table_data = []
    all_ring_ids = sorted(all_ring_ids)
    
    # 预测结构列表（排除参考结构）
    predict_structures = [s for s in results.keys() if s != reference_name]
    
    for ring_id in all_ring_ids:
        row = [ring_id]
        ref_grid = reference_rings.get(ring_id, ['NA', 'NA', 'NA'])
        row.append(ref_grid)
        # 对每个预测结构添加 grid unit 和与参考结构的距离
        for struct in predict_structures:
            pred_grid = results[struct].get(ring_id, ['NA', 'NA', 'NA'])
            if 'NA' in ref_grid or 'NA' in pred_grid:
                distance = 'NA'
            else:
                distance = calculate_distance(ref_grid, pred_grid)
            row.append(pred_grid)
            row.append(distance)
        table_data.append(row)
    
    # 构造 DataFrame 列名
    columns = ["Ring Identifier", "Reference Grid Unit"]
    for struct in predict_structures:
        columns.append(f"{struct} Grid Unit")
        columns.append(f"{struct} Distance")
    
    df = pd.DataFrame(table_data, columns=columns)
    output_csv = os.path.join(cache_folder, "ligand_ring_grid_unit_comparison_with_distance.csv")
    df.to_csv(output_csv, index=False)
    print(f"环的 grid unit 比较结果已保存到 {output_csv}")
    
    return results

# 如需将此功能注册为 PyMOL 命令，请取消下一行注释：
# cmd.extend("find_ligand_ring_grid_units", find_ligand_ring_grid_units)
