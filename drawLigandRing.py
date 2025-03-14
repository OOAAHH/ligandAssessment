# -*- coding: utf-8 -*-
# output csv
# Author: HaoSun
# Date: 2025-03-12 19:08
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
draw_ring_centers_gridunits.py

本插件对给定选择（sele）中的非氢原子计算环的质量中心，
将质心转换为 grid unit 坐标（公式： n = floor((c + spacing/2)/spacing) ），
并以伪原子的形式绘制环的质心，同时为每个质心在偏移位置标注带有前缀 "MC" 的 grid unit 坐标标签。

新增参数：
  atom_size    : 控制伪原子球体的大小，默认值 0.5

用法:
    draw_ring_centers_gridunits sele, [name], [spacing], [label_offset], [atom_size]
    参数:
      sele         : PyMOL 中的选择条件（例如 "ligand"）
      name         : 生成的伪原子对象名称，默认 "ringCenters"
      spacing      : 网格单元间距，默认 2.0 Å
      label_offset : 标签偏移值（沿 Z 方向），默认 0.5 Å
      atom_size    : 伪原子球体的大小，默认 0.5
"""

import math
from pymol import cmd
import networkx as nx

def draw_ring_centers_gridunits(sele, name="ringCenters", spacing=2.0, label_offset=0.5, atom_size=0.5):
    spacing = float(spacing)
    label_offset = float(label_offset)
    atom_size = float(atom_size)
    model = cmd.get_model(sele)
    if not model.atom:
        print("未找到原子，请检查选择条件。")
        return

    # 定义常见元素的原子质量和共价半径
    atomic_masses = {
        'C': 12.01,
        'N': 14.01,
        'O': 16.00,
        'P': 30.97,
        'S': 32.07,
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

    # 构建非氢原子字典，key 为 atom.index
    atoms_dict = {}
    for atom in model.atom:
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

    if len(atoms_dict) < 3:
        print("不足3个非氢原子，无法构成环。")
        return

    # 构建图：节点为 atom.index
    G = nx.Graph()
    for idx in atoms_dict.keys():
        G.add_node(idx)

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
            r_i = atoms_dict[idx_i]['radius']
            r_j = atoms_dict[idx_j]['radius']
            if dist < (r_i + r_j + tolerance):
                G.add_edge(idx_i, idx_j)

    # 利用 networkx 的 cycle_basis 检测所有简单环
    cycles = nx.cycle_basis(G)
    if not cycles:
        print("未检测到环。")
        return

    # grid unit 坐标转换函数：对于 spacing=2，
    # 网格单元 n 对应真实坐标区间 [2*n - 1, 2*n + 1)
    def convert_to_grid_unit(coordinate, spacing):
        return [math.floor((c + spacing/2) / spacing) for c in coordinate]

    # 清除旧的对象
    cmd.delete(name)
    cmd.delete(name + "_label")

    # 对每个环计算质量中心，转换为 grid unit 坐标，并创建伪原子和标签
    for cycle in cycles:
        if len(cycle) < 3:
            continue
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
        grid_unit = convert_to_grid_unit(center_of_mass, spacing)
        # 在标签前增加 "MC" 前缀
        label_str = f"MC({grid_unit[0]},{grid_unit[1]},{grid_unit[2]})"
        
        # 在质心位置创建伪原子（用于显示球体，并通过 vdw 参数调整大小）
        cmd.pseudoatom(object=name, pos=center_of_mass, label="", vdw=atom_size)
        # 在偏移位置创建伪原子，用于显示标签，避免被球体遮挡
        offset_pos = [center_of_mass[0], center_of_mass[1], center_of_mass[2] + label_offset]
        cmd.pseudoatom(object=name + "_label", pos=offset_pos, label=label_str)

    # 显示球体和标签
    cmd.show("spheres", name)
    cmd.show("labels", name + "_label")
    
    # 设置标签文字的颜色为明亮黄色，加粗字体
    cmd.set("label_color", "yellow", name + "_label")
    cmd.set("label_font_id", 7, name + "_label")
    
    print("已绘制环的质量中心，并标记带有 'MC' 前缀的 grid unit 坐标（标签已偏移且为黄色加粗）。")

# 注册为 PyMOL 命令
cmd.extend("draw_ring_centers_gridunits", draw_ring_centers_gridunits)
