# -*- coding: utf-8 -*-
# output csv
# Author: HaoSun
# Date: 2025-02-12 14:55
import json
import csv
import argparse

def format_contact(contact):
    """
    根据 contact 字典生成格式化字符串：
    "model1,atom1-model2,atom2; (x1,y1,z1) - (x2,y2,z2); calculated_distance Å"
    坐标和距离均保留三位小数。
    """
    model1 = contact["model1"]
    atom1 = contact["atom1"]
    coord1 = contact["coord1"]
    model2 = contact["model2"]
    atom2 = contact["atom2"]
    coord2 = contact["coord2"]
    distance = contact["calculated_distance"]
    coord1_str = f"({coord1[0]:.3f},{coord1[1]:.3f},{coord1[2]:.3f})"
    coord2_str = f"({coord2[0]:.3f},{coord2[1]:.3f},{coord2[2]:.3f})"
    return f"{model1},{atom1}-{model2},{atom2}; {coord1_str} - {coord2_str}; {distance:.3f} Å"

def compare_contacts(xray_contact, predict_contact, tol=0.1):
    """
    比较 x‑ray 与预测 contact 的 calculated_distance，
    若两者之差小于 tol，则返回 "✔️ 完全匹配"，否则返回 "❌ 部分偏移"。
    """
    diff = abs(xray_contact["calculated_distance"] - predict_contact["calculated_distance"])
    return "✔️ 完全匹配" if diff < tol else "❌ 部分偏移"

def process_json_to_table(json_file, xray_identifier):
    """
    从 JSON 文件中读取多个组数据，每个组格式为：
      {
         "pdb_file": "6GZR.1-10_09.pdb",
         "contacts": [ { contact1 }, { contact2 }, ... ]
      }
    根据外部参数 xray_identifier 判断哪个组为 x‑ray 组，其它组均为预测组。
    
    对于每个预测组，生成多行数据，格式为：
      [contact序号, predict pdb name, x-ray, predict, 匹配状态]
    其中：
      - contact序号：如果索引在 x‑ray 联系中存在，则为 "contact{i+1}"；否则为 "extra{j}"（j 从1开始）。
      - x‑ray 列：如果对应 x‑ray 组中有该 contact，则格式化显示，否则为空。
      - predict 列：显示该预测组中该 contact 的格式化信息（若存在）。
      - 匹配状态：若两边均存在，则比较 calculated_distance，否则为空。
      
    返回：所有行数据的列表。
    """
    with open(json_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    # 分离 x‑ray 组和预测组
    xray_group = None
    predict_groups = []
    for group in data:
        pdb_file = group.get("pdb_file", "")
        if xray_identifier in pdb_file:
            xray_group = group
        else:
            predict_groups.append(group)
    
    if xray_group is None:
        raise ValueError(f"未找到包含 '{xray_identifier}' 的 x‑ray 组")
    if len(predict_groups) == 0:
        raise ValueError("未找到预测组")
    
    xray_contacts = xray_group.get("contacts", [])
    num_xray = len(xray_contacts)
    
    # 对于每个预测组，生成行数据
    rows = []
    for pred_group in predict_groups:
        pred_pdb = pred_group.get("pdb_file", "")
        pred_contacts = pred_group.get("contacts", [])
        num_pred = len(pred_contacts)
        total = max(num_xray, num_pred)
        for i in range(total):
            # 如果 i < num_xray，则 contact_id = "contact{i+1}"；否则 "extra{j}"
            if i < num_xray:
                contact_id = f"contact {i+1}"
                xray_str = format_contact(xray_contacts[i])
            else:
                contact_id = f"extra {i - num_xray + 1}"
                xray_str = ""
            # 预测数据
            if i < num_pred:
                pred_str = format_contact(pred_contacts[i])
            else:
                pred_str = ""
            # 匹配状态：如果两边都有数据则比较，否则为空
            if i < num_xray and i < num_pred:
                match_status = compare_contacts(xray_contacts[i], pred_contacts[i])
            else:
                match_status = ""
            row = [contact_id, pred_pdb, xray_str, pred_str, match_status]
            rows.append(row)
    
    # 表头
    header = ["contact序号", "predict pdb name", "x-ray", "predict", "匹配状态"]
    return header, rows

def write_csv(header, rows, output_csv):
    """
    将表格数据写入 CSV 文件。
    """
    with open(output_csv, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"CSV 文件已保存至 {output_csv}")

def convert_json_to_csv(json_file, config_file, output_csv):
    """
    将 JSON 转换为 CSV 格式。
    """
    # 从 config 文件读取 select 信息
    with open(config_file, 'r') as f:
        lines = f.readlines()
        xray_identifier = lines[0].strip()  # 获取 x-ray 组的标识符

    header, rows = process_json_to_table(json_file, xray_identifier)
    write_csv(header, rows, output_csv)
