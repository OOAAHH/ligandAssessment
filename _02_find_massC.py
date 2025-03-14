# -*- coding: utf-8 -*-
# 对于非氢原子, 程序会先尝试在字典 atomic_masses 中查找该原子的质量.     #
# 如果没有找到对应的质量(即原子符号不在 atomic_masses 字典中),          #
# 那么就会使用默认值 12.0 作为该原子的质量. 这种方式允许你在不明确指定     #
# 所有可能原子类型的情况下仍能进行计算, 但默认值可能需要根据实际情况调整    #
# output mass center json
# Author: HaoSun
# Date: 2025-03-06 17:16
import os
import json
from pymol import cmd

def read_config(config_file):
    """
    读取配置文件, 返回参考结构名称和ligand选择条件. 
    
    配置文件要求：
      第一行：参考结构名称(不带扩展名)
      第二行：ligand选择条件(例如 "chain A" 或 "resn LIG")
    """
    with open(config_file, 'r') as f:
        lines = f.readlines()
    if len(lines) < 2:
        raise ValueError("配置文件至少需要两行：第一行参考结构名称, 第二行ligand选择条件")
    reference_name = lines[0].strip()
    select_ligand_chain = lines[1].strip()
    return reference_name, select_ligand_chain

def load_structures(input_folder):
    """
    加载指定文件夹中所有PDB文件的完整路径, 并返回排序后的列表. 
    """
    structures = []
    for root, _, files in os.walk(input_folder):
        for file in files:
            if file.lower().endswith(".pdb"):
                structures.append(os.path.join(root, file))
    return sorted(structures)

def get_reference_structure(structures, reference_name):
    """
    在结构列表中查找参考结构文件, 通过比较文件基本名称(不含扩展名)匹配. 
    """
    for structure in structures:
        basename = os.path.basename(structure).replace(".pdb", "")
        if basename == reference_name:
            return structure
    raise ValueError(f"未在输入文件夹中找到参考结构：{reference_name}")

def compute_center_of_mass(selection):
    """
    计算指定选择中重原子(非氢原子)的质心(加权中心). 
    
    使用常见原子的质量进行加权计算, 返回 [x, y, z] 三维坐标列表, 
    若选择中无重原子则返回 None. 
    """
    model = cmd.get_model(selection)
    if not model.atom:
        return None

    # 定义常见原子的原子质量(单位：amu)
    atomic_masses = {
        'C': 12.01,
        'N': 14.01,
        'O': 16.00,
        'P': 30.97,
        'S': 32.07,
        # 可根据需要添加更多元素
    }
    total_mass = 0.0
    weighted_sum = [0.0, 0.0, 0.0]
    
    for atom in model.atom:
        # 尝试获取原子元素符号, 若无则使用 elem 属性
        element = getattr(atom, 'symbol', None)
        if element is None:
            element = atom.elem
        # 排除氢原子
        if element.upper() == 'H':
            continue
        mass = atomic_masses.get(element.upper(), 12.0)  # 默认值为12.0
        total_mass += mass
        weighted_sum[0] += atom.coord[0] * mass
        weighted_sum[1] += atom.coord[1] * mass
        weighted_sum[2] += atom.coord[2] * mass

    if total_mass == 0:
        return None
    center_of_mass = [weighted_sum[0] / total_mass,
                      weighted_sum[1] / total_mass,
                      weighted_sum[2] / total_mass]
    return center_of_mass

def find_ligand_mass_center(input_folder, config_file):
    """
    对所有 predict 结构的 ligand 区域计算重原子质心. 
    
    使用与之前代码相同的 ligand 选择逻辑(依据配置文件第二行的条件), 
    并跳过参考结构(配置文件第一行指定). 
    
    返回一个字典, key 为 predict 结构的基本名称, value 为 ligand 重原子质心 [x, y, z]. 
    同时将结果保存到 JSON 文件 "ligand_mass_center.json" 中. 
    """
    # 读取配置文件内容
    reference_name, select_ligand_chain = read_config(config_file)

    # 检查 cache 文件夹
    cache_folder = os.path.join(input_folder, "cache")
    if not os.path.exists(cache_folder):
        print("Run super_rnas() first")

    # 加载所有结构
    structures = load_structures(cache_folder)

    # 加载参考结构
    reference_file = get_reference_structure(structures, reference_name)
    print(f"加载参考结构: {reference_name}")
    cmd.load(reference_file, reference_name)
    
    # 用于存储每个 predict 结构 ligand 重原子质心的结果
    results = {}
    
    # 遍历所有结构, 跳过参考结构
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
        
        # 使用与 align 逻辑相同的 ligand 选择条件
        cmd.select("ligand", f"{basename} and {select_ligand_chain}")
        
        # 计算 ligand 区域中重原子的质心
        mass_center = compute_center_of_mass("ligand")
        if mass_center is not None:
            print(f"{basename} 的 ligand 重原子质心: {mass_center}")
            results[basename] = mass_center
        else:
            print(f"{basename} 的 ligand 区域中未找到重原子")
            results[basename] = None
        
        # 清理已加载的结构和选择
        cmd.delete(basename)
        cmd.delete("ligand")
    
    # 将结果保存到 JSON 文件中
    output_file = os.path.join(cache_folder, "ligand_mass_center.json")
    with open(output_file, "w") as f:
        json.dump(results, f, indent=4)
    print(f"ligand 重原子质心结果已保存到 {output_file}")
    
    return results


def find_massC(input_folder, config_file):
    """
    主函数：执行ligand的对齐, 并输出RMSD值到JSON文件. 
    
    参数:
      input_folder : 包含PDB文件的文件夹路径
      config_file  : 包含参考结构名称和ligand选择条件的配置文件路径
    """
    find_ligand_mass_center(input_folder, config_file)