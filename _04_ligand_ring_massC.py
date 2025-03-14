# -*- coding: utf-8 -*-
# 计算 ligand 中每个非氢原子环的重原子质心
# Author: HaoSun
# Date: 2025-03-12
# input super的结构和reference，output ligand_ring_mass_center.json
#
# 说明：
# 1. 对每个结构，根据配置文件中指定的 ligand 选择条件提取非氢原子；
# 2. 根据原子之间的距离（利用各元素的共价半径和一定容差）构建连接图；
# 3. 利用 networkx.cycle_basis 算法检测环；
# 4. 对每个环计算质量加权的质心，结果保存到 JSON 文件中。

import os
import json
import math
from pymol import cmd
import networkx as nx

def read_config(config_file):
    """
    读取配置文件，返回参考结构名称和 ligand 选择条件。
    
    配置文件要求：
      第一行：参考结构名称（不带扩展名，可用于区分参考和预测结构，可选）
      第二行：ligand 选择条件（例如 "chain A" 或 "resn LIG"）
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
    加载指定文件夹中所有 PDB 文件的完整路径，并返回排序后的列表。
    """
    structures = []
    for root, _, files in os.walk(input_folder):
        for file in files:
            if file.lower().endswith(".pdb"):
                structures.append(os.path.join(root, file))
    return sorted(structures)

def compute_ring_mass_centers(selection):
    """
    对指定选择（通常为 ligand 区域）中的非氢原子，
    1. 构建原子间的连接图（根据距离判断是否存在化学键）；
    2. 利用 networkx.cycle_basis 检测图中的所有简单环；
    3. 对每个环计算质量加权的质心，返回一个列表，列表中每个元素为字典：
         { "atoms": [atom_name1, atom_name2, ...],
           "atom_indices": [index1, index2, ...],
           "center_of_mass": [x, y, z] }
    若选择中无非氢原子或无环，则返回空列表。
    """
    model = cmd.get_model(selection)
    if not model.atom:
        return []

    # 定义常见元素的原子质量（单位：amu）和共价半径（单位：Å）
    atomic_masses = {
        'C': 12.01,
        'N': 14.01,
        'O': 16.00,
        'P': 30.97,
        'S': 32.07,
        # 可根据需要补充
    }
    covalent_radii = {
        'C': 0.76,
        'N': 0.71,
        'O': 0.66,
        'P': 1.07,
        'S': 1.05,
    }
    default_mass = 12.0
    default_radius = 0.77
    tolerance = 0.4  # 距离容差

    # 构建一个字典保存非氢原子信息，key 为 atom.index，value 为字典（name, element, coord, mass, radius）
    atoms_dict = {}
    for atom in model.atom:
        # 尝试获取原子元素符号
        element = getattr(atom, 'symbol', None)
        if element is None:
            element = atom.elem
        if element.upper() == 'H':
            continue
        mass = atomic_masses.get(element.upper(), default_mass)
        radius = covalent_radii.get(element.upper(), default_radius)
        atoms_dict[atom.index] = {
            'name': atom.name,
            'element': element.upper(),
            'coord': atom.coord,
            'mass': mass,
            'radius': radius
        }

    # 若不足3个原子，不可能形成环
    if len(atoms_dict) < 3:
        return []

    # 构建图：节点为 atom.index
    G = nx.Graph()
    for idx in atoms_dict.keys():
        G.add_node(idx)

    # 对于每对原子，判断两原子间距离是否小于 (radius1 + radius2 + tolerance)
    atom_indices = list(atoms_dict.keys())
    n = len(atom_indices)
    for i in range(n):
        for j in range(i+1, n):
            idx_i = atom_indices[i]
            idx_j = atom_indices[j]
            coord_i = atoms_dict[idx_i]['coord']
            coord_j = atoms_dict[idx_j]['coord']
            dx = coord_i[0] - coord_j[0]
            dy = coord_i[1] - coord_j[1]
            dz = coord_i[2] - coord_j[2]
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            # 判断是否“成键”
            r_i = atoms_dict[idx_i]['radius']
            r_j = atoms_dict[idx_j]['radius']
            if dist < (r_i + r_j + tolerance):
                G.add_edge(idx_i, idx_j)

    # 利用 networkx 的 cycle_basis 寻找所有简单环
    cycles = nx.cycle_basis(G)
    ring_centers = []
    for cycle in cycles:
        if len(cycle) < 3:
            continue  # 环至少需要3个原子
        total_mass = 0.0
        weighted_sum = [0.0, 0.0, 0.0]
        atom_names = []
        for idx in cycle:
            atom_info = atoms_dict[idx]
            m = atom_info['mass']
            total_mass += m
            weighted_sum[0] += atom_info['coord'][0] * m
            weighted_sum[1] += atom_info['coord'][1] * m
            weighted_sum[2] += atom_info['coord'][2] * m
            atom_names.append(atom_info['name'])
        if total_mass == 0:
            continue
        center_of_mass = [weighted_sum[0] / total_mass,
                          weighted_sum[1] / total_mass,
                          weighted_sum[2] / total_mass]
        ring_centers.append({
            "atoms": atom_names,
            "atom_indices": cycle,
            "center_of_mass": center_of_mass
        })

    return ring_centers

def find_ligand_ring_mass_center(input_folder, config_file):
    """
    对 cache 文件夹中的所有结构，根据配置文件中的 ligand 选择条件，
    计算 ligand 区域中每个非氢原子环的重原子质心。
    
    处理过程：
      1. 读取配置文件，获取参考结构名称（可用于排除）和 ligand 选择条件；
      2. 加载 cache 文件夹中所有 PDB 文件，依次载入结构，选择 ligand 区域；
      3. 调用 compute_ring_mass_centers 计算每个结构中 ligand 内所有环的质心；
      4. 将每个结构的结果保存到字典中，并写入 JSON 文件 "ligand_ring_mass_center.json"。
    """
    # 读取配置文件
    reference_name, select_ligand_chain = read_config(config_file)

    # 检查 cache 文件夹
    cache_folder = os.path.join(input_folder, "cache")
    if not os.path.exists(cache_folder):
        print("Cache 文件夹不存在，请先运行相关计算生成 cache 文件夹。")
        return {}

    # 加载所有结构（可以选择跳过参考结构，如果需要）
    structures = load_structures(cache_folder)
    
    results = {}
    for structure in structures:
        basename = os.path.basename(structure).replace(".pdb", "")
        # 如需跳过参考结构，可取消下面注释
        # if basename == reference_name:
        #     print(f"跳过参考结构: {reference_name}")
        #     continue

        print(f"处理结构 {basename} ...")
        # 加载结构至 PyMOL，并使用配置中指定的 ligand 选择条件
        cmd.load(structure, basename)
        cmd.select("ligand", f"{basename} and {select_ligand_chain}")
        
        # 计算 ligand 内各环的重原子质心
        ring_centers = compute_ring_mass_centers("ligand")
        if ring_centers:
            print(f"{basename} 中检测到 {len(ring_centers)} 个环。")
            results[basename] = ring_centers
        else:
            print(f"{basename} 中未检测到非氢原子环。")
            results[basename] = []
        
        # 清理加载的结构和选择
        cmd.delete(basename)
        cmd.delete("ligand")
    
    # 保存结果到 JSON 文件
    output_file = os.path.join(cache_folder, "ligand_ring_mass_center.json")
    with open(output_file, "w") as f:
        json.dump(results, f, indent=4)
    print(f"Ligand 每个环的重原子质心结果已保存到 {output_file}")
    
    return results


# cmd.extend("find_ligand_ring_mass_center", find_ligand_ring_mass_center)
