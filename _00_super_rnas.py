# -*- coding: utf-8 -*-
# Author: HaoSun
# Date: 2025-01-27 16:30
import os
from pymol import cmd

def read_config(config_file):
    """
    读取配置文件，返回参考结构名称和RNA选择条件。
    
    配置文件要求：
      第一行：参考结构名称（不带扩展名）
      第三行：RNA选择条件（例如 "chain A"）
    """
    with open(config_file, 'r') as f:
        lines = f.readlines()
    if len(lines) < 3:
        raise ValueError("配置文件至少需要三行：第一行参考结构名称，第三行RNA选择条件")
    reference_name = lines[0].strip()
    select_rna_chain = lines[2].strip()
    return reference_name, select_rna_chain

def load_structures(input_folder):
    """
    加载指定文件夹中所有PDB文件的完整路径，并返回排序后的列表。
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

def load_and_align_structures(input_folder, config_file):
    """
    加载所有结构，加载参考结构并对齐其他结构到参考结构。
    对齐后的结构保存到 input_folder/cache 文件夹中。
    
    参数:
      input_folder : 包含PDB文件的文件夹路径
      config_file  : 包含参考结构名称和RNA选择条件的配置文件路径
    """
    # 读取配置文件
    reference_name, select_rna_chain = read_config(config_file)
    
    # 加载所有结构
    structures = load_structures(input_folder)
    
    # 查找参考结构文件
    reference_file = get_reference_structure(structures, reference_name)
    
    # 创建缓存文件夹
    cache_folder = os.path.join(input_folder, "cache")
    if not os.path.exists(cache_folder):
        os.makedirs(cache_folder)
    
    # 加载参考结构，指定对象名称为 reference_name
    print(f"加载参考结构: {reference_name} ...")
    cmd.load(reference_file, reference_name)
    # 保存参考结构到缓存文件夹
    ref_save_path = os.path.join(cache_folder, f"{reference_name}.pdb")
    cmd.save(ref_save_path, reference_name)
    print(f"保存参考结构到缓存文件夹: {reference_name} ...")
    
    # 遍历其他结构并对齐
    for structure in structures:
        basename = os.path.basename(structure).replace(".pdb", "")
        if basename == reference_name:
            print(f"跳过参考结构: {reference_name}")
            continue
        
        print(f"对齐结构 {basename} 到参考结构 {reference_name} ...")
        # 加载当前结构，并指定对象名称为 basename
        cmd.load(structure, basename)
        
        # 通过选择创建RNA区域
        cmd.select("rna1", f"{basename} and {select_rna_chain}")
        cmd.select("rna2", f"{reference_name} and {select_rna_chain}")
        
        # 使用 super 进行结构对齐
        cmd.super("rna1", "rna2")
        
        # 保存对齐后的结构到缓存文件夹
        save_path = os.path.join(cache_folder, f"{basename}_super.pdb")
        cmd.save(save_path, basename)
        
        # 删除当前加载的结构和选择，保持工作空间整洁
        cmd.delete(basename)
        cmd.delete("rna1")
        cmd.delete("rna2")
        
        print(f"完成对齐结构: {basename}")

def super_rnas(input_folder, config_file):
    """
    主函数：执行所有结构的加载和对齐操作。
    
    参数:
      input_folder : 包含PDB文件的文件夹路径
      config_file  : 包含参考结构名称和RNA选择条件的配置文件路径
    """
    load_and_align_structures(input_folder, config_file)
