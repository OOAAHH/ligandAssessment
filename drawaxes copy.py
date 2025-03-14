# -*- coding: utf-8 -*-
# output csv
# Author: HaoSun
# Date: 2025-03-12 16:47
"""
drawaxes_gridunits.py

本插件根据所有加载的模型计算包围盒，并按照新的逻辑：
1. 遍历所有原子，分别获取X、Y、Z方向的最小值和最大值；
2. 根据每个方向上真实坐标与网格单元的关系（对于 spacing=2 时，网格单元 n 对应真实坐标区间 [2*n-1, 2*n+1)），
   计算出每个方向的网格下标范围；
3. 计算每个方向上的小立方体真实边界（必为奇数位置，当 spacing=2）；
4. 在每个方向上，仅在这些边界处绘制直线作为坐标轴线（包括正负方向）。
同时，为非氢原子添加 grid unit 坐标标签。

用法 drawaxes_gridunits sele
"""

import math
from pymol.cgo import BEGIN, LINES, VERTEX, END, COLOR
from pymol import cmd

# 更新后的坐标转换函数
def convert_to_grid_unit(coordinate, spacing=2.0):
    """
    对于 spacing=2 的情况，满足：
        网格单元 n 对应真实区间 [2*n-1, 2*n+1)
    因此，给定坐标 c，可用公式： n = floor((c + spacing/2) / spacing)
    """
    grid_unit = []
    for c in coordinate:
        grid_unit.append(math.floor((c + spacing/2) / spacing))
    return grid_unit

def drawaxes_gridunits(sele, name="axesGrid", spacing=2.0):
    spacing = float(spacing)

    all_atoms = cmd.get_model("all").atom
    if not all_atoms:
        print("没有找到任何原子。")
        return

    # 提取所有原子的 X, Y, Z 坐标
    xs = [atom.coord[0] for atom in all_atoms]
    ys = [atom.coord[1] for atom in all_atoms]
    zs = [atom.coord[2] for atom in all_atoms]

    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)

    # 根据公式：n = floor((c + spacing/2)/spacing)
    # 计算每个方向上最小和最大的网格下标（注意最大值用 ceil 后减1保证包含边界）
    n_min_x = math.floor((min_x + spacing/2) / spacing)
    n_max_x = math.ceil((max_x + spacing/2) / spacing) - 1

    n_min_y = math.floor((min_y + spacing/2) / spacing)
    n_max_y = math.ceil((max_y + spacing/2) / spacing) - 1

    n_min_z = math.floor((min_z + spacing/2) / spacing)
    n_max_z = math.ceil((max_z + spacing/2) / spacing) - 1

    # 根据网格下标计算真实坐标边界：
    # 对于网格单元 n，真实边界为：[n*spacing - spacing/2, n*spacing + spacing/2]
    aligned_min_x = n_min_x * spacing - spacing/2
    aligned_max_x = n_max_x * spacing + spacing/2
    aligned_min_y = n_min_y * spacing - spacing/2
    aligned_max_y = n_max_y * spacing + spacing/2
    aligned_min_z = n_min_z * spacing - spacing/2
    aligned_max_z = n_max_z * spacing + spacing/2

    # 构造每个方向上的边界列表（间隔为 spacing，且对于 spacing=2 时边界均为奇数）
    x_boundaries = []
    cur = aligned_min_x
    while cur <= aligned_max_x + 1e-6:
        x_boundaries.append(cur)
        cur += spacing

    y_boundaries = []
    cur = aligned_min_y
    while cur <= aligned_max_y + 1e-6:
        y_boundaries.append(cur)
        cur += spacing

    z_boundaries = []
    cur = aligned_min_z
    while cur <= aligned_max_z + 1e-6:
        z_boundaries.append(cur)
        cur += spacing

    # 开始构造 CGO 以绘制坐标轴线（直线仅绘制在各方向的边界上，即奇数位置）
    cgo_grid = [BEGIN, LINES]

    # 绘制平行于 X 轴的直线：
    # 固定 Y 和 Z 坐标（均取边界值），X 坐标从 X 方向边界的最小值绘制到最大值
    for y in y_boundaries:
        for z in z_boundaries:
            cgo_grid.extend([COLOR, 0.8, 0.5, 0.5])  # 红色：表示 X 轴方向的线
            cgo_grid.extend([VERTEX, aligned_min_x, y, z])
            cgo_grid.extend([VERTEX, aligned_max_x, y, z])

    # 绘制平行于 Y 轴的直线：
    # 固定 X 和 Z 坐标，Y 坐标从 Y 方向边界的最小值绘制到最大值
    for x in x_boundaries:
        for z in z_boundaries:
            cgo_grid.extend([COLOR, 0.5, 0.8, 0.5])  # 绿色：表示 Y 轴方向的线
            cgo_grid.extend([VERTEX, x, aligned_min_y, z])
            cgo_grid.extend([VERTEX, x, aligned_max_y, z])

    # 绘制平行于 Z 轴的直线：
    # 固定 X 和 Y 坐标，Z 坐标从 Z 方向边界的最小值绘制到最大值
    for x in x_boundaries:
        for y in y_boundaries:
            cgo_grid.extend([COLOR, 0.5, 0.5, 0.8])  # 蓝色：表示 Z 轴方向的线
            cgo_grid.extend([VERTEX, x, y, aligned_min_z])
            cgo_grid.extend([VERTEX, x, y, aligned_max_z])

    cgo_grid.append(END)
    cmd.load_cgo(cgo_grid, f"{name}_grid")

    # 根据新的转换函数为非氢原子添加 grid unit 坐标标签
    cmd.select("ligand_heavy", f"{sele} and not elem H")
    heavy_model = cmd.get_model("ligand_heavy")
    for atom in heavy_model.atom:
        grid_unit = convert_to_grid_unit(atom.coord, spacing)
        label_str = f"({grid_unit[0]},{grid_unit[1]},{grid_unit[2]})"
        sel_atom = f"ligand_heavy and index {atom.index}"
        cmd.alter(sel_atom, f'label="{label_str}"')
    cmd.label("ligand_heavy", "label")
    print("已为 ligand heavy atoms 添加 grid unit 坐标标签。")

cmd.extend("drawaxes_gridunits", drawaxes_gridunits)
