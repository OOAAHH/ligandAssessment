# -*- coding: utf-8 -*-
"""
draw_grid_path.py

本插件要求在 PyMOL 中选中恰好两个原子，
从第一个原子出发，沿着 x、y、z 三个方向构造 Manhattan 路径，
路径采用虚线显示（dash/gap 模式）且使用圆柱体绘制以增强醒目效果；
在每段路径的中点（沿指定方向偏移2Å）处添加加粗橙色 label 显示跨越的 grid 数量，
并在 label 中增加下划线效果；
同时，在起点和终点处添加显示 grid unit 坐标的 label（偏移3Å，显示为白色），
在绘制的坐标系正中间输出路径 grid unit 总和（三个方向的和，显示为 cyan）。
此外，根据两个原子的位置，计算并绘制正好覆盖它们的最小 grid 坐标轴。

为防止绘制多个 path 时名称冲突，增加了一个自定义名称参数。

用法示例：
    select mysel, (index 10+index 50)
    draw_grid_path sele, name="massCentre", spacing=2.0, line_width=0.6
"""

import math
from pymol.cgo import CYLINDER, BEGIN, LINES, VERTEX, END, COLOR, LINEWIDTH
from pymol import cmd

def sanitize_name(name):
    """仅保留字母和数字，其它字符替换为下划线，并转换为小写"""
    return "".join([c if c.isalnum() else "_" for c in name]).lower()

def generate_dashed_segments(p_start, p_end, dash_length, gap_length):
    """
    给定起始点和终点，按照 dash_length 与 gap_length 生成虚线段列表。
    返回格式为 [(s1, e1), (s2, e2), ...]。
    """
    segments = []
    dx = p_end[0] - p_start[0]
    dy = p_end[1] - p_start[1]
    dz = p_end[2] - p_start[2]
    seg_length = math.sqrt(dx*dx + dy*dy + dz*dz)
    if seg_length == 0:
        return segments
    ux, uy, uz = dx/seg_length, dy/seg_length, dz/seg_length
    t = 0.0
    while t < seg_length:
        dash_start = t
        dash_end = min(t + dash_length, seg_length)
        s = (p_start[0] + ux*dash_start, p_start[1] + uy*dash_start, p_start[2] + uz*dash_start)
        e = (p_start[0] + ux*dash_end, p_start[1] + uy*dash_end, p_start[2] + uz*dash_end)
        segments.append((s, e))
        t += dash_length + gap_length
    return segments

def convert_to_grid_unit(coord, spacing):
    """
    将一个坐标转换为 grid unit 坐标，公式为：
      grid_unit = floor((c + spacing/2) / spacing)
    返回一个 (grid_x, grid_y, grid_z) 的元组。
    """
    return (math.floor((coord[0] + spacing/2) / spacing),
            math.floor((coord[1] + spacing/2) / spacing),
            math.floor((coord[2] + spacing/2) / spacing))

def add_cylinder_line(cgo_list, p_start, p_end, radius, color):
    """
    将一段 dash 段绘制为圆柱体，颜色 color 为 (r, g, b)
    """
    cgo_list.extend([CYLINDER,
                     p_start[0], p_start[1], p_start[2],
                     p_end[0], p_end[1], p_end[2],
                     radius,
                     color[0], color[1], color[2],
                     color[0], color[1], color[2]])

def draw_grid_path(sele, name="grid_path", spacing=2.0, line_width=5.0, dash_length=None, gap_length=None):
    """
    从选中两个原子构造 Manhattan 路径（依次改变 x、y、z 坐标）：
      - 路径虚线采用 dash_length 与 gap_length 分段绘制（默认均为 spacing/2），
        每段以圆柱体显示，radius 可根据 line_width 调整；
      - 每段中点沿固定方向偏移2Å，并以加粗橙色 label 显示跨越的 grid 数量，
        且在 label 中增加下划线效果；
      - 在起点和终点处额外添加显示 grid unit 坐标的 label（偏移3Å，显示为白色）；
      - 在绘制的 grid 坐标轴正中间输出路径 grid unit 总和（三个方向的和，显示为 cyan）；
      - 同时，根据两个原子确定的最小包围盒绘制 grid 坐标轴（采用原有颜色方案）。
    
    参数：
      sele       : 包含两个原子的选择名
      name       : 自定义名称前缀，防止多个 path 时对象覆盖（将自动清洗），默认 "grid_path"
      spacing    : grid 间距，默认2.0Å
      line_width : 用于计算 dash 绘制时的“粗度”，默认5.0
      dash_length: 虚线每段长度（默认 spacing/2）
      gap_length : 虚线间隙长度（默认 spacing/2）
    """
    # 对名称进行清洗，确保生成合法对象名
    name = sanitize_name(name)

    # 确保数值类型
    spacing = float(spacing)
    line_width = float(line_width)
    if dash_length is None:
        dash_length = spacing / 2.0
    else:
        dash_length = float(dash_length)
    if gap_length is None:
        gap_length = spacing / 2.0
    else:
        gap_length = float(gap_length)
        
    # 为 dash 段圆柱体设置半径（此处设为 line_width 的 0.2 倍，可根据需要调整）
    cylinder_radius = line_width * 0.2

    model = cmd.get_model(sele)
    if len(model.atom) != 2:
        print("请确保选择中恰好有两个原子。")
        return

    # 获取两个原子的坐标
    atom1 = model.atom[0]
    atom2 = model.atom[1]
    x1, y1, z1 = atom1.coord
    x2, y2, z2 = atom2.coord

    # 计算 Manhattan 路径节点：p1 -> (x2, y1, z1) -> (x2, y2, z1) -> p2
    p1 = (x1, y1, z1)
    p_mid1 = (x2, y1, z1)
    p_mid2 = (x2, y2, z1)
    p2 = (x2, y2, z2)

    # 计算各段跨越的 grid 数量（四舍五入取整）
    grid_count_x = int(round(abs(x2 - x1) / spacing))
    grid_count_y = int(round(abs(y2 - y1) / spacing))
    grid_count_z = int(round(abs(z2 - z1) / spacing))
    grid_total = grid_count_x + grid_count_y + grid_count_z

    # 构造 dash 圆柱体的 CGO 对象，逐段绘制
    dashed_cgo = []
    dash_color = (1.0, 0.5, 0.0)  # 橙色
    def add_dashed_line(p_start, p_end):
        segs = generate_dashed_segments(p_start, p_end, dash_length, gap_length)
        for seg in segs:
            s, e = seg
            add_cylinder_line(dashed_cgo, s, e, cylinder_radius, dash_color)
    add_dashed_line(p1, p_mid1)
    add_dashed_line(p_mid1, p_mid2)
    add_dashed_line(p_mid2, p2)
    cmd.load_cgo(dashed_cgo, f"{name}_path")

    # 根据两个原子计算用于绘制 grid 坐标轴的最小包围盒
    min_x, max_x = min(x1, x2), max(x1, x2)
    min_y, max_y = min(y1, y2), max(y1, y2)
    min_z, max_z = min(z1, z2), max(z1, z2)

    n_min_x = math.floor((min_x + spacing/2) / spacing)
    n_max_x = math.ceil((max_x + spacing/2) / spacing) - 1
    n_min_y = math.floor((min_y + spacing/2) / spacing)
    n_max_y = math.ceil((max_y + spacing/2) / spacing) - 1
    n_min_z = math.floor((min_z + spacing/2) / spacing)
    n_max_z = math.ceil((max_z + spacing/2) / spacing) - 1

    aligned_min_x = n_min_x * spacing - spacing/2
    aligned_max_x = n_max_x * spacing + spacing/2
    aligned_min_y = n_min_y * spacing - spacing/2
    aligned_max_y = n_max_y * spacing + spacing/2
    aligned_min_z = n_min_z * spacing - spacing/2
    aligned_max_z = n_max_z * spacing + spacing/2

    def build_boundaries(aligned_min, aligned_max, spacing):
        boundaries = []
        cur = aligned_min
        while cur <= aligned_max + 1e-6:
            boundaries.append(cur)
            cur += spacing
        return boundaries

    x_boundaries = build_boundaries(aligned_min_x, aligned_max_x, spacing)
    y_boundaries = build_boundaries(aligned_min_y, aligned_max_y, spacing)
    z_boundaries = build_boundaries(aligned_min_z, aligned_max_z, spacing)

    # 构造 grid 坐标轴 CGO（采用原有颜色方案）
    grid_cgo = [BEGIN, LINES]
    # 平行 X 轴的线（固定 y 和 z）
    for y in y_boundaries:
        for z in z_boundaries:
            grid_cgo.extend([COLOR, 0.8, 0.5, 0.5,
                             VERTEX, aligned_min_x, y, z,
                             VERTEX, aligned_max_x, y, z])
    # 平行 Y 轴的线（固定 x 和 z）
    for x in x_boundaries:
        for z in z_boundaries:
            grid_cgo.extend([COLOR, 0.5, 0.8, 0.5,
                             VERTEX, x, aligned_min_y, z,
                             VERTEX, x, aligned_max_y, z])
    # 平行 Z 轴的线（固定 x 和 y）
    for x in x_boundaries:
        for y in y_boundaries:
            grid_cgo.extend([COLOR, 0.5, 0.5, 0.8,
                             VERTEX, x, y, aligned_min_z,
                             VERTEX, x, y, aligned_max_z])
    grid_cgo.append(END)
    cmd.load_cgo(grid_cgo, f"{name}_axes")

    # 计算各段 Manhattan 路径的中点（用于路径 label），偏移仍为2Å
    mid_x = [(p1[0] + p_mid1[0]) / 2.0, (p1[1] + p_mid1[1]) / 2.0, (p1[2] + p_mid1[2]) / 2.0]
    mid_y = [(p_mid1[0] + p_mid2[0]) / 2.0, (p_mid1[1] + p_mid2[1]) / 2.0, (p_mid1[2] + p_mid2[2]) / 2.0]
    mid_z = [(p_mid2[0] + p2[0]) / 2.0, (p_mid2[1] + p2[1]) / 2.0, (p_mid2[2] + p2[2]) / 2.0]
    # 固定偏移方向：X段 label 沿正 Y 偏移2Å, Y段 label 沿正 Z 偏移2Å, Z段 label 沿正 X 偏移2Å
    offset1 = (0, 2.0, 0)
    offset2 = (0, 0, 2.0)
    offset3 = (2.0, 0, 0)
    mid1_offset = [mid_x[i] + offset1[i] for i in range(3)]
    mid2_offset = [mid_y[i] + offset2[i] for i in range(3)]
    mid3_offset = [mid_z[i] + offset3[i] for i in range(3)]

    # 为各段生成 label 文本，并在文本后增加下划线（下划线长度与字符数相同）
    label1_text = f"X: {grid_count_x}"
    label2_text = f"Y: {grid_count_y}"
    label3_text = f"Z: {grid_count_z}"
    underline1 = "_" * len(label1_text)
    underline2 = "_" * len(label2_text)
    underline3 = "_" * len(label3_text)
    label1 = f"{label1_text}\n{underline1}"
    label2 = f"{label2_text}\n{underline2}"
    label3 = f"{label3_text}\n{underline3}"

    # 创建描述各段路径长度的 label（用 pseudoatom 实现），局部设置其颜色为橙色和加粗
    cmd.pseudoatom(object=f"{name}_label1", pos=mid1_offset, label=label1)
    cmd.pseudoatom(object=f"{name}_label2", pos=mid2_offset, label=label2)
    cmd.pseudoatom(object=f"{name}_label3", pos=mid3_offset, label=label3)
    cmd.set("label_color", "orange", f"obj {name}_label1 or obj {name}_label2 or obj {name}_label3")
    cmd.set("label_font_id", 7, f"obj {name}_label1 or obj {name}_label2 or obj {name}_label3")

    # 对起点和终点添加 grid unit 坐标 label，偏移改为3Å（此处默认沿正 Y 方向偏移3Å）
    start_grid = convert_to_grid_unit(p1, spacing)
    end_grid = convert_to_grid_unit(p2, spacing)
    label_start = f"Start: {start_grid}"
    label_end = f"End: {end_grid}"
    start_offset = (0, 3.0, 0)
    end_offset = (0, 3.0, 0)
    start_pos = [p1[i] + start_offset[i] for i in range(3)]
    end_pos = [p2[i] + end_offset[i] for i in range(3)]
    cmd.pseudoatom(object=f"{name}_start", pos=start_pos, label=label_start)
    cmd.pseudoatom(object=f"{name}_end", pos=end_pos, label=label_end)
    cmd.set("label_color", "white", f"obj {name}_start or obj {name}_end")

    # 计算 grid 坐标轴中心位置，并在此处显示三个方向 grid 数量的总和（summary）
    center_x = (aligned_min_x + aligned_max_x) / 2.0
    center_y = (aligned_min_y + aligned_max_y) / 2.0
    center_z = (aligned_min_z + aligned_max_z) / 2.0
    summary_pos = (center_x, center_y, center_z)
    summary_label = f"Summary: {grid_total}"
    cmd.pseudoatom(object=f"{name}_summary", pos=summary_pos, label=summary_label)
    cmd.set("label_color", "cyan", f"obj {name}_summary")
    cmd.set("label_font_id", 7, f"obj {name}_summary")

    print(f"路径绘制成功：X方向 {grid_count_x} 格，Y方向 {grid_count_y} 格，Z方向 {grid_count_z} 格，总计 {grid_total} 格。")

# 将命令添加到 PyMOL 中
cmd.extend("draw_grid_path", draw_grid_path)
