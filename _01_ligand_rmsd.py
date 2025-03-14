# -*- coding: utf-8 -*-
# output rmsd json
# Author: HaoSun
# Date: 2025-02-03 21:10
import os
import json
from pymol import cmd

def read_config(config_file):
    """
    读取配置文件，返回参考结构名称和ligand选择条件。
    
    配置文件要求：
      第一行：参考结构名称（不带扩展名）
      第二行：ligand选择条件（例如 "chain A" 或 "resn LIG"）
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
    加载所有结构，对齐非参考结构的ligand部分到参考结构的ligand部分，
    并使用rms_cur方法计算RMSD值，最后将结果输出到JSON文件中。
    
    每一个predict结构对应到x‑ray结构的RMSD值，
    结果以predict结构文件的基本名称为key，RMSD值为value进行分组。
    
    参数:
      input_folder : 包含PDB文件的文件夹路径
      config_file  : 包含参考结构名称和ligand选择条件的配置文件路径
    """
    # 读取配置文件内容
    reference_name, select_ligand_chain = read_config(config_file)
    
    # 检查或创建cache文件夹
    cache_folder = os.path.join(input_folder, "cache")
    if not os.path.exists(cache_folder):
        print("Run super_rnas() to create cache folder: {cache_folder}")
    
    # 加载所有PDB结构
    structures = load_structures(cache_folder)
    
    # 查找并加载参考结构
    reference_file = get_reference_structure(structures, reference_name)
    print(f"加载参考结构: {reference_name}")
    cmd.load(reference_file, reference_name)
    ref_save_path = os.path.join(cache_folder, f"{reference_name}.pdb")
    cmd.save(ref_save_path, reference_name)
    
    # 用于保存结果的字典，key为predict结构的基本名称，value为RMSD值
    results = {}
    
    # 遍历其他结构进行对齐
    for structure in structures:
        basename = os.path.basename(structure).replace(".pdb", "")
        if basename == reference_name:
            print(f"跳过参考结构: {reference_name}")
            continue
        
        print(f"对齐结构 {basename} 到参考结构 {reference_name} ...")
        # 加载当前结构
        cmd.load(structure, basename)
        
        # 选择ligand区域
        cmd.select("ligand1", f"{basename} and {select_ligand_chain}")
        cmd.select("ligand2", f"{reference_name} and {select_ligand_chain}")
        
        # 使用rms_cur方法计算RMSD值
        rmsd_value = cmd.rms_cur("ligand1", "ligand2")
        if rmsd_value is not None:
            print(f"使用rms_cur方法计算 {basename} 的RMSD: {rmsd_value:.2f}")
            results[basename] = rmsd_value
        else:
            print(f"对齐 {basename} 未能获得RMSD值")
            results[basename] = None
        
        # 删除当前加载的结构和选择，保持PyMOL工作空间整洁
        cmd.delete(basename)
        cmd.delete("ligand1")
        cmd.delete("ligand2")
        
        print(f"完成对齐结构: {basename}")
    
    # 将结果输出到JSON文件中
    output_file = os.path.join(input_folder, "rmsd_results.json")
    with open(output_file, "w") as f:
        json.dump(results, f, indent=4)
    print(f"结果已保存到 {output_file}")

def ligand_rmsd(input_folder, config_file):
    """
    主函数：执行ligand的对齐，并输出RMSD值到JSON文件。
    
    参数:
      input_folder : 包含PDB文件的文件夹路径
      config_file  : 包含参考结构名称和ligand选择条件的配置文件路径
    """
    load_and_align_structures(input_folder, config_file)
