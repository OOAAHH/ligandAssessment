# -*- coding: utf-8 -*-
import numpy as np
from pymol.cgo import BEGIN, LINES, VERTEX, END, COLOR
from pymol import cmd

def findmincuboid_pi(sele, name="mincuboid_pi"):
    """
    计算给定选择 (sele) 中所有非氢原子（重原子）的最小体积定向包围盒，
    要求该包围盒有两个面平行于ligand最佳π平面（即π平面的法向量为包围盒的一个坐标轴）。
    
    方法：
      1. 从选择中提取所有非氢原子坐标；
      2. 计算质心并对坐标进行PCA，得到三个主轴，其中最小特征值对应的特征向量
         作为最佳π平面法向量；
      3. 构造旋转矩阵：令第三轴（z轴）等于该法向量，另外两个轴取自PCA结果并用叉乘确保右手系，
         这样转换后xy平面就与ligand的π平面平行；
      4. 将点转换到新坐标系下，计算各轴的min/max构造包围盒；
      5. 将包围盒顶点转换回原始坐标系，并利用CGO以红色线框显示。
      
    用法示例：
      PyMOL命令行中输入：
         findmincuboid_pi sele, mincuboid_pi_obj
      或在Python中调用：
         cmd.findmincuboid_pi("sele", name="mincuboid_pi_obj")
    """
    # 1. 提取非氢原子坐标
    model = cmd.get_model(sele)
    coords = []
    for atom in model.atom:
        element = getattr(atom, 'symbol', None)
        if element is None:
            element = atom.elem
        if element.upper() == 'H':
            continue
        coords.append(atom.coord)
    if not coords:
        print("选择区域中未找到重原子。")
        return
    points = np.array(coords)
    
    # 2. 计算质心及PCA
    center = points.mean(axis=0)
    points_centered = points - center
    cov = np.cov(points_centered, rowvar=False)
    eigvals, eigvecs = np.linalg.eigh(cov)
    # np.linalg.eigh返回的是升序排列的特征值，最小的对应最佳π平面法向量
    # 为了方便构造坐标系，我们按降序排列
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    # 设：
    #   v1: 第一主轴（最大方差方向）
    #   v2: 第二主轴
    #   v3: 第三主轴（最小方差方向，作为最佳π平面法向量）
    v1 = eigvecs[:, 0]
    v3 = eigvecs[:, 2]  # π平面的法向量
    # 令v2' = v3 x v1，确保右手系
    v2 = np.cross(v3, v1)
    v2 = v2 / np.linalg.norm(v2)
    # 重新正交化：更新v1 = v2 x v3
    v1 = np.cross(v2, v3)
    v1 = v1 / np.linalg.norm(v1)
    
    # 3. 构造旋转矩阵 R，其列向量为 v1, v2, v3
    R = np.column_stack((v1, v2, v3))
    
    # 4. 坐标转换
    rotated = np.dot(points_centered, R)
    min_rot = rotated.min(axis=0)
    max_rot = rotated.max(axis=0)
    
    # 构造新坐标系下包围盒的8个顶点（利用min_rot和max_rot）
    # 其中 v0 = [min_x, min_y, min_z], v1 = [max_x, min_y, min_z], etc.
    v0 = np.array([min_rot[0], min_rot[1], min_rot[2]])
    v1_ = np.array([max_rot[0], min_rot[1], min_rot[2]])
    v2_ = np.array([max_rot[0], max_rot[1], min_rot[2]])
    v3 = np.array([min_rot[0], max_rot[1], min_rot[2]])
    v4 = np.array([min_rot[0], min_rot[1], max_rot[2]])
    v5 = np.array([max_rot[0], min_rot[1], max_rot[2]])
    v6 = np.array([max_rot[0], max_rot[1], max_rot[2]])
    v7 = np.array([min_rot[0], max_rot[1], max_rot[2]])
    
    # 将包围盒顶点从新坐标系转换回原始坐标系： x = (vertex in rotated) * R^T + center
    V0 = np.dot(v0, R.T) + center
    V1 = np.dot(v1_, R.T) + center
    V2 = np.dot(v2_, R.T) + center
    V3 = np.dot(v3, R.T) + center
    V4 = np.dot(v4, R.T) + center
    V5 = np.dot(v5, R.T) + center
    V6 = np.dot(v6, R.T) + center
    V7 = np.dot(v7, R.T) + center
    
    # 5. 定义长方体12条边，并构造CGO对象绘制
    edges = [
        (V0, V1), (V1, V2), (V2, V3), (V3, V0),  # 底面
        (V4, V5), (V5, V6), (V6, V7), (V7, V4),  # 顶面
        (V0, V4), (V1, V5), (V2, V6), (V3, V7)   # 侧面
    ]
    
    cgo = [BEGIN, LINES, COLOR, 1.0, 0.0, 0.0]
    for edge in edges:
        for vertex in edge:
            cgo.extend([VERTEX, float(vertex[0]), float(vertex[1]), float(vertex[2])])
    cgo.append(END)
    
    cmd.load_cgo(cgo, name)
    print("定向包围盒对象 '%s' 已加载。" % name)

# 将该函数注册为PyMOL命令
cmd.extend("findmincuboid_pi", findmincuboid_pi)
