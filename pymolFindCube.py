# -*- coding: utf-8 -*-
import numpy as np
import os
from pymol.cgo import BEGIN, LINES, VERTEX, END, COLOR
from pymol import cmd
import time

def findmincuboid_pi(sele, name="mincuboid_pi", view_face=False, ray=False, save_png=False, png_filename=None):
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
      6. 可选: 调整视角面向立方体的面，渲染场景并保存为PNG图像。
      
    参数:
      sele: 用于计算包围盒的选择对象
      name: 生成的包围盒对象名称
      view_face: 是否调整视角面向立方体的最大面
      ray: 是否进行光线追踪渲染
      save_png: 是否保存PNG图像
      png_filename: PNG文件名(可选)，若未指定则使用name作为文件名
      
    用法示例：
      PyMOL命令行中输入：
         findmincuboid_pi sele, mincuboid_pi_obj
         findmincuboid_pi sele, mincuboid_pi_obj, view_face=1, ray=1, save_png=1
      或在Python中调用：
         cmd.findmincuboid_pi("sele", name="mincuboid_pi_obj")
         cmd.findmincuboid_pi("sele", name="mincuboid_pi_obj", view_face=True, ray=True, save_png=True, png_filename="my_view.png")
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
    
    # 6. 根据需要调整视角、渲染和保存图像
    if view_face:
        # 立方体CGO对象不能用于orient，直接对原始选择进行orient
        cmd.orient(sele, animate=0)
        
        # 确保立方体也在视图中
        cmd.zoom(sele, buffer=2.0)
        
        print(f"已将视角调整为面向选中结构，立方体可见。")
        
        # 对场景进行光线追踪渲染
        if ray:
            cmd.ray()
            print("已完成光线追踪渲染。")
        
        # 保存为PNG图像
        if save_png:
            if png_filename is None:
                png_filename = f"{name}.png"
            cmd.png(png_filename)
            print(f"已将图像保存为: {os.path.abspath(png_filename)}")

# 将该函数注册为PyMOL命令
cmd.extend("findmincuboid_pi", findmincuboid_pi)

# 添加一个新函数，通过伪原子创建可以使用orient命令的立方体对象
def findmincuboid_pi_with_atoms(sele, name="mincuboid_pi_atoms", view_face=False, ray=False, save_png=False, png_filename=None):
    """
    与findmincuboid_pi功能相同，但使用伪原子创建立方体，使其可以使用orient命令。
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
    
    # 清除之前的对象（如果存在）
    cmd.delete(name)
    
    # 创建一个新的空分子对象
    cmd.create(name, None)
    
    # 使用伪原子创建立方体顶点
    vertices = [V0, V1, V2, V3, V4, V5, V6, V7]
    for i, vertex in enumerate(vertices):
        # 将NumPy数组转换为元组以避免格式问题
        x, y, z = float(vertex[0]), float(vertex[1]), float(vertex[2])
        cmd.pseudoatom(
            object=name,
            selection='',
            name=f'V{i}',
            resn='BOX',
            resi=str(i),
            chain='X',
            segi='CUBE',
            elem='C',
            vdw=0.5,
            hetatm=1,
            b=0.0,
            q=1.0,
            color='red',
            label='',
            pos=[x, y, z]  # 使用Python列表而不是NumPy数组
        )
    
    # 定义立方体的边（键）
    edges = [
        (0, 1), (1, 2), (2, 3), (3, 0),  # 底面
        (4, 5), (5, 6), (6, 7), (7, 4),  # 顶面
        (0, 4), (1, 5), (2, 6), (3, 7)   # 侧面连接
    ]
    
    # 使用bond命令创建键而不是distance对象
    for v1, v2 in edges:
        cmd.bond(f"{name} and resi {v1} and name V{v1}", f"{name} and resi {v2} and name V{v2}")
    
    # 设置显示风格
    cmd.show("spheres", name)
    cmd.show("sticks", name)
    cmd.set("sphere_scale", 0.3, name)
    cmd.set("stick_radius", 0.1, name)
    cmd.color("red", name)
    
    print(f"定向包围盒对象 '{name}' 已加载(使用伪原子表示)。")
    
    # 根据需要调整视角、渲染和保存图像
    if view_face:
        # 使用更可靠的方式调整视角
        # 首先确保对象存在
        if cmd.count_atoms(name) > 0:
            # 使用完整的对象名称进行orient
            cmd.orient(f"({name})", animate=0)
            print(f"已将视角调整为面向立方体 '{name}'。")
        else:
            print(f"警告: 立方体对象 '{name}' 不存在或为空。")
            # 尝试对原始选择进行orient
            cmd.orient(sele, animate=0)
            print(f"已将视角调整为面向原始选择 '{sele}'。")
        
        # 添加短暂延迟，确保旋转操作完全生效
        time.sleep(0.5)
        
        # 对场景进行光线追踪渲染
        if ray:
            cmd.ray()
            print("已完成光线追踪渲染。")
        
        # 保存为PNG图像
        if save_png:
            if png_filename is None:
                png_filename = f"{name}.png"
            cmd.png(png_filename)
            print(f"已将图像保存为: {os.path.abspath(png_filename)}")

# 将新函数注册为PyMOL命令
cmd.extend("findmincuboid_pi_with_atoms", findmincuboid_pi_with_atoms)

# 添加一个新函数，专门用于调整立方体视角
def orient_cube(cube_name, animate=0):
    """
    调整视角面向指定的立方体对象。
    
    参数:
      cube_name: 立方体对象的名称
      animate: 是否使用动画效果（默认为0，表示不使用动画）
    """
    # 检查对象是否存在
    if cmd.count_atoms(cube_name) > 0:
        # 使用完整的对象名称进行orient
        cmd.orient(f"({cube_name})", animate=animate)
        print(f"已将视角调整为面向立方体 '{cube_name}'。")

    else:
        print(f"错误: 立方体对象 '{cube_name}' 不存在或为空。")

# 将新函数注册为PyMOL命令
cmd.extend("orient_cube", orient_cube)
