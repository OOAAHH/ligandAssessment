# -*- coding: utf-8 -*-
# output ligand cube
# Author: HaoSun
# Date: 2025-01-27 19:51
# Edit: 2025-03-07 14:28
# Notes: 更改寻找cube的逻辑，现在是基于最佳pai平面的算法。
# 新增ligand cube的可视化，在pymol中显示ligand cube：ligand_workflow/pymolFindCube.py

import os
import json
import numpy as np
from pymol import cmd

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

def compute_oriented_cube_from_heavy_atoms(selection):
    """
    计算指定选择中重原子（非氢原子）的定向最小包围立方体，
    要求立方体的两个面与ligand最佳π平面平行。
    
    算法步骤：
      1. 提取选择中所有重原子坐标；
      2. 计算点云的质心，并对点云做PCA；取方差最小的特征向量作为π平面法向量，
         同时选取最大方向以及通过叉乘获得正交方向构造旋转矩阵 R，
         使得 R 的第三列即为π平面法向量，前两列构成平面内的正交基；
      3. 将所有点平移至以质心为原点后，转换到新坐标系 rotated = (points - center) * R，
         在新坐标系中计算各轴最小、最大值，再取最大跨度作为cube边长；
      4. 以新坐标系下包围盒中心构造立方体顶点（中心 ± half_side），
         然后转换回原始坐标系： vertex = center + (vertex_rotated)*R^T；
      5. 返回包含8个顶点坐标的列表，每个顶点为 [x, y, z]。
    """
    model = cmd.get_model(selection)
    if not model.atom:
        return None

    coords = []
    for atom in model.atom:
        element = getattr(atom, 'symbol', None)
        if element is None:
            element = atom.elem
        if element.upper() == 'H':
            continue
        coords.append(atom.coord)
    
    if not coords:
        return None
    
    points = np.array(coords)
    # 计算质心
    center = points.mean(axis=0)
    points_centered = points - center
    
    # PCA计算：协方差矩阵和特征分解（特征值按升序排列，最小值对应π平面法向量）
    cov = np.cov(points_centered, rowvar=False)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)  # 升序排列
    # 取最小特征值对应的特征向量作为π平面法向量
    v_pi = eigvecs[:, order[0]]
    # 取最大特征值对应的方向作为平面内一方向
    v_long = eigvecs[:, order[-1]]
    # 通过叉乘获得另一正交方向
    v_ortho = np.cross(v_pi, v_long)
    # 归一化
    v_pi /= np.linalg.norm(v_pi)
    v_long /= np.linalg.norm(v_long)
    v_ortho /= np.linalg.norm(v_ortho)
    
    # 构造旋转矩阵 R，其列向量为 v_long, v_ortho, v_pi
    # 这样新坐标系的 z 轴与ligand π平面法向量一致，新坐标系的 xy 平面与π平面平行
    R = np.column_stack((v_long, v_ortho, v_pi))
    
    # 将点转换到新坐标系
    rotated = np.dot(points_centered, R)
    
    # 计算在新坐标系下的各轴最小、最大值
    min_rot = rotated.min(axis=0)
    max_rot = rotated.max(axis=0)
    spans = max_rot - min_rot
    # 取最大跨度作为cube的边长
    side = np.max(spans)
    half_side = side / 2.0
    # 计算新坐标系下cube中心
    center_rot = (min_rot + max_rot) / 2.0
    
    # 计算cube 8个顶点：中心 ± half_side 各轴组合
    corners = []
    for dx in (-half_side, half_side):
        for dy in (-half_side, half_side):
            for dz in (-half_side, half_side):
                vertex_rot = center_rot + np.array([dx, dy, dz])
                # 转换回原始坐标系
                vertex = center + np.dot(vertex_rot, R.T)
                corners.append(vertex.tolist())
    
    return corners

def find_ligand_cube(input_folder, config_file):
    """
    对所有 predict 结构的 ligand 区域，利用最佳π平面思路计算定向最小包围cube，
    返回每个结构对应的8个顶点坐标。
    
    使用配置文件中第二行的 ligand 选择条件，跳过参考结构（配置文件第一行指定）。
    
    返回字典：key 为 predict 结构基本名称，value 为该结构 ligand 的cube 8个顶点坐标列表，
    并保存结果到 JSON 文件 "ligand_cube.json"（保存在 cache 文件夹中）。
    """
    # 读取配置文件内容
    reference_name, select_ligand_chain = read_config(config_file)

    # 检查或创建 cache 文件夹（如有其他用途）
    cache_folder = os.path.join(input_folder, "cache")
    if not os.path.exists(cache_folder):
        print("Run super_rnas() first")
    # 加载所有结构（这里假设已在 cache 文件夹中存放predict结构）
    structures = load_structures(cache_folder)

    # 加载参考结构（用于后续排除）
    reference_file = get_reference_structure(structures, reference_name)
    print(f"加载参考结构: {reference_name}")
    cmd.load(reference_file, reference_name)
    
    # 用于存储每个 predict 结构 ligand 包围cube 8个顶点的结果
    results = {}
    
    # 遍历所有结构，跳过参考结构
    for structure in structures:
        basename = os.path.basename(structure).replace(".pdb", "")
        if basename == reference_name:
            print(f"跳过参考结构: {reference_name}")
            continue
        
        print(f"处理结构 {basename} ...")
        # 加载 predict 结构
        cmd.load(structure, basename)
        
        # 选择 ligand 区域（使用配置文件中指定的选择条件）
        cmd.select("ligand", f"{basename} and {select_ligand_chain}")
        
        # 计算 ligand 重原子的定向cube顶点（利用最佳π平面思路）
        cube_corners = compute_oriented_cube_from_heavy_atoms("ligand")
        if cube_corners is not None:
            print(f"{basename} 的 ligand cube 顶点:")
            for pt in cube_corners:
                print(f"  {pt}")
            results[basename] = cube_corners
        else:
            print(f"{basename} 的 ligand 区域中未找到重原子")
            results[basename] = None
        
        # 清理加载的结构和选择
        cmd.delete(basename)
        cmd.delete("ligand")
    
    # 将结果保存到 JSON 文件中（保存在 cache 文件夹下）
    output_file = os.path.join(cache_folder, "ligand_cube.json")
    with open(output_file, "w") as f:
        json.dump(results, f, indent=4)
    print(f"Ligand cube 顶点结果已保存到 {output_file}")
    
    return results

# 重命名函数，便于调用
def compute_oriented_cube_from_heavy_atoms(selection):
    """
    利用最佳π平面思路，计算给定选择中重原子的定向最小包围cube，
    要求cube的两个面与ligand最佳π平面平行。
    """
    return compute_oriented_cube_from_heavy_atoms_inner(selection)

def compute_oriented_cube_from_heavy_atoms_inner(selection):
    # （参考上面的函数实现）
    model = cmd.get_model(selection)
    if not model.atom:
        return None

    coords = []
    for atom in model.atom:
        element = getattr(atom, 'symbol', None)
        if element is None:
            element = atom.elem
        if element.upper() == 'H':
            continue
        coords.append(atom.coord)
    
    if not coords:
        return None
    
    points = np.array(coords)
    center = points.mean(axis=0)
    points_centered = points - center

    cov = np.cov(points_centered, rowvar=False)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)
    # 最小特征值对应的特征向量作为最佳π平面法向量
    v_pi = eigvecs[:, order[0]]
    # 最大特征值对应的方向
    v_long = eigvecs[:, order[-1]]
    v_ortho = np.cross(v_pi, v_long)
    v_pi /= np.linalg.norm(v_pi)
    v_long /= np.linalg.norm(v_long)
    v_ortho /= np.linalg.norm(v_ortho)
    
    # 构造旋转矩阵 R：使得新坐标系的z轴为 v_pi，xy平面与π平面平行
    R = np.column_stack((v_long, v_ortho, v_pi))
    
    rotated = np.dot(points_centered, R)
    min_rot = rotated.min(axis=0)
    max_rot = rotated.max(axis=0)
    spans = max_rot - min_rot
    side = np.max(spans)
    half_side = side / 2.0
    center_rot = (min_rot + max_rot) / 2.0

    corners = []
    for dx in (-half_side, half_side):
        for dy in (-half_side, half_side):
            for dz in (-half_side, half_side):
                vertex_rot = center_rot + np.array([dx, dy, dz])
                vertex = center + np.dot(vertex_rot, R.T)
                corners.append(vertex.tolist())
    return corners

def ligand_cube(input_folder, config_file):
    """
    主函数：执行所有 predict 结构中 ligand 的定向cube计算，
    输出结果保存到JSON文件中。
    
    参数:
      input_folder : 包含PDB文件的文件夹路径
      config_file  : 包含参考结构名称和ligand选择条件的配置文件路径
    """
    find_ligand_cube(input_folder, config_file)
