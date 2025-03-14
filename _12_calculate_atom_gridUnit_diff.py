# -*- coding: utf-8 -*-
# output grid unit distance
# 读取重原子的真实坐标的json，输出grid unit坐标的
# Author: HaoSun
# Date: 2025-03-12 11:09 (更新于新版逻辑)
import os
import json
import numpy as np
import pandas as pd
from pymol import cmd
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
    读取配置文件，返回参考结构名称和ligand选择条件。
    """
    with open(config_file, 'r') as f:
        lines = f.readlines()
    if len(lines) < 2:
        raise ValueError("配置文件至少需要两行：第一行参考结构名称，第二行ligand选择条件")
    reference_name = lines[0].strip()
    select_ligand_chain = lines[1].strip()
    return reference_name, select_ligand_chain

def load_structures(input_folder):
    """
    加载指定文件夹中的所有PDB文件，并返回排序后的列表。
    """
    structures = []
    for root, _, files in os.walk(input_folder):
        for file in files:
            if file.lower().endswith(".pdb"):
                structures.append(os.path.join(root, file))
    return sorted(structures)

def get_reference_structure(structures, reference_name):
    """
    在结构列表中查找参考结构文件，通过比较文件基本名称（不含扩展名）匹配。
    """
    for structure in structures:
        basename = os.path.basename(structure).replace(".pdb", "")
        if basename == reference_name:
            return structure
    raise ValueError(f"未在输入文件夹中找到参考结构：{reference_name}")

def calculate_distance(grid_unit_ref, grid_unit_predict):
    """
    计算参考结构和预测结构的 grid unit 坐标的距离
    计算公式：|x_ref - x_predict| + |y_ref - y_predict| + |z_ref - z_predict|
    """
    distance = sum(abs(np.array(grid_unit_ref) - np.array(grid_unit_predict)))
    return distance

def find_ligand_grid_units(input_folder, config_file):
    """
    对所有 predict 结构的 ligand 区域，利用 convert_to_grid_unit(coordinate)
    返回每个结构对应的重原子 grid unit 坐标，并计算与参考结构的距离。

    使用配置文件中第二行的 ligand 选择条件，跳过参考结构（配置文件第一行指定）。

    返回的 grid unit 坐标以字典形式返回，键为结构名称，值为重原子 grid unit 坐标列表，同时保存为 JSON 文件。
    """
    # 读取配置文件内容
    reference_name, select_ligand_chain = read_config(config_file)

    # 检查 cache 文件夹
    cache_folder = os.path.join(input_folder, "cache")
    if not os.path.exists(cache_folder):
        print("Run super_rnas() first")

    # 加载所有结构
    structures = load_structures(cache_folder)

    # 加载参考结构（用于后续排除）
    reference_file = get_reference_structure(structures, reference_name)
    print(f"加载参考结构: {reference_name}")
    cmd.load(reference_file, reference_name)

    # 用于存储每个 predict 结构 ligand 的 grid unit 坐标
    results = {}
    reference_grid_units = {}

    # 计算参考结构的 grid unit 坐标
    cmd.select("ligand_ref", f"{reference_name} and {select_ligand_chain}")
    ligand_ref_coords = cmd.get_model("ligand_ref")
    ligand_ref_coordinates = {}
    for atom in ligand_ref_coords.atom:
        element = getattr(atom, 'symbol', None)
        if element is None:
            element = atom.elem
        if element.upper() == 'H':
            continue
        ligand_ref_coordinates[atom.name] = atom.coord

    for atom_name, coord in ligand_ref_coordinates.items():
        grid_unit = convert_to_grid_unit(coord)  # 使用新逻辑计算 grid unit 坐标
        reference_grid_units[atom_name] = grid_unit

    # 遍历所有结构，跳过参考结构
    for structure in structures:
        basename = os.path.basename(structure).replace(".pdb", "")
        '''
        if basename == reference_name:
            print(f"跳过参考结构: {reference_name}")
            continue
        '''
        
        print(f"处理结构 {basename} ...")
        # 加载 predict 结构
        cmd.load(structure, basename)
        
        # 选择 ligand 区域（使用配置文件中指定的选择条件）
        cmd.select("ligand", f"{basename} and {select_ligand_chain}")
        
        # 获取 ligand 重原子坐标，排除氢原子
        ligand_coords = cmd.get_model("ligand")
        ligand_coordinates = {}
        for atom in ligand_coords.atom:
            element = getattr(atom, 'symbol', None)
            if element is None:
                element = atom.elem
            if element.upper() == 'H':
                continue
            ligand_coordinates[atom.name] = atom.coord
        
        # 计算 ligand 重原子 grid unit 坐标
        ligand_grid_units = {}
        for atom_name, coord in ligand_coordinates.items():
            grid_unit = convert_to_grid_unit(coord)  # 使用新逻辑计算 grid unit 坐标
            ligand_grid_units[atom_name] = grid_unit
        
        # 将转换后的结果存储到字典中
        results[basename] = ligand_grid_units
        
        # 清理加载的结构和选择
        cmd.delete(basename)
        cmd.delete("ligand")
    
    # 将所有结构的 ligand grid unit 结果保存到 JSON 文件中（保存在 cache 文件夹下）
    output_file = os.path.join(cache_folder, "ligand_massAtom_gridUnit.json")
    with open(output_file, "w") as f:
        json.dump(results, f, indent=4)
    print(f"Ligand grid unit结果已保存到 {output_file}")

    # 创建一个表格来比较每个结构的重原子 grid unit 结果与参考结构的差异
    table_data = []
    atom_names = list(reference_grid_units.keys())  # 使用参考结构中的原子名称作为横轴
    for atom_name in atom_names:
        row = [atom_name]  # 重原子名称
        row.append(reference_grid_units.get(atom_name, ['NA', 'NA', 'NA']))  # 参考结构的 grid unit
        for structure in structures:
            basename = os.path.basename(structure).replace(".pdb", "")
            if basename == reference_name:
                continue
            grid_unit_predict = results.get(basename, {}).get(atom_name, ['NA', 'NA', 'NA'])
            distance = calculate_distance(reference_grid_units[atom_name], grid_unit_predict)
            row.append(grid_unit_predict)  # 预测结构的 grid unit 坐标
            row.append(distance)  # 与参考结构的距离
        table_data.append(row)

    # 构造 DataFrame 并保存为 CSV 文件
    df = pd.DataFrame(
        table_data,
        columns=["Atom Name", "Reference Grid Unit"] +
                [f"{os.path.basename(structure).replace('.pdb', '')} Grid Unit" for structure in structures if os.path.basename(structure).replace(".pdb", "") != reference_name] +
                [f"{os.path.basename(structure).replace('.pdb', '')} Distance" for structure in structures if os.path.basename(structure).replace(".pdb", "") != reference_name]
    )
    df.to_csv(os.path.join(cache_folder, "ligand_grid_unit_comparison_with_distance.csv"), index=False)

    return results