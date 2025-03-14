# -*- coding: utf-8 -*-
# output contacts json
# Author: HaoSun
# Date: 2025-02-12 14:22
import os
import json
import pymol
from pymol import cmd
import math

def polarpairs(sel1, sel2, cutoff=4.0, angle=54.0, name='', state=1, quiet=1):
    cutoff = float(cutoff)
    quiet = int(quiet)
    state = int(state)
    if angle == 'default':
        angle = cmd.get('h_bond_max_angle', cmd.get_object_list(sel1)[0])
    angle = float(angle)
    mode = 1 if angle > 0 else 0
    x = cmd.find_pairs('(%s) and donors' % sel1, '(%s) and acceptors' % sel2,
            state, state,
            cutoff=cutoff, mode=mode, angle=angle) + \
        cmd.find_pairs('(%s) and acceptors' % sel1, '(%s) and donors' % sel2,
            state, state,
            cutoff=cutoff, mode=mode, angle=angle)
    x = sorted(set(x))
    if not quiet:
        print('Settings: cutoff=%.1f Å, angle=%.1f°' % (cutoff, angle))
        print('Found %d polar contacts' % (len(x)))
    if len(name) > 0:
        for p in x:
            cmd.distance(name, '(%s%s)' % p[0], '(%s%s)' % p[1])
    return x

def get_atom_coordinates(model_id, atom_index):
    selection = f"{model_id} and index {atom_index}"
    model = cmd.get_model(selection)
    if len(model.atom) == 0:
        raise ValueError(f"No atom found with index {atom_index} in model {model_id}")
    atom_coord = model.atom[0].coord
    return tuple(round(coord, 3) for coord in atom_coord)

def calculate_distance(coord1, coord2):
    distance = math.sqrt((coord1[0] - coord2[0])**2 +
                         (coord1[1] - coord2[1])**2 +
                         (coord1[2] - coord2[2])**2)
    return round(distance, 3)

def process_contacts_and_calculate_distance(contact_data):
    results = []
    for contact in contact_data:
        model1, atom1 = contact[0]
        model2, atom2 = contact[1]
        original_distance = contact[2]
        
        coord1 = get_atom_coordinates(model1, atom1)
        coord2 = get_atom_coordinates(model2, atom2)
        
        calc_distance = calculate_distance(coord1, coord2)
        
        contact_result = {
            'model1': model1,
            'atom1': atom1,
            'coord1': coord1,
            'model2': model2,
            'atom2': atom2,
            'coord2': coord2,
            'calculated_distance': calc_distance,
            'original_distance': original_distance
        }
        results.append(contact_result)

    return results

def process_directory(input_folder, config_file, output_json):
    """
    读取config文件，解析并批量处理指定文件夹中的所有PDB文件，计算接触信息。
    """
    # 读取config文件中的选择信息
    with open(config_file, 'r') as f:
        lines = f.readlines()
        select_ligand = lines[1].strip()  # 获取ligand选择条件
        select_rna_chain = lines[2].strip()  # 获取RNA链选择条件

    all_results = []

    for pdb_file in os.listdir(input_folder):
        if pdb_file.endswith('.pdb'):
            pdb_path = os.path.join(input_folder, pdb_file)
            print(f"Processing {pdb_file}...")

            # 加载PDB文件
            cmd.load(pdb_path)

            cmd.select("ligand", select_ligand)
            cmd.select("rna_chain", select_rna_chain)

            contacts = polarpairs('rna_chain', 'ligand', angle=63, cutoff=3.5)

            contact_data = []
            for p in contacts:
                orig_dist = cmd.get_distance('(%s and index %d)' % (p[0][0], p[0][1]),
                                            '(%s and index %d)' % (p[1][0], p[1][1]))
                orig_dist = round(orig_dist, 3)
                contact_data.append((p[0], p[1], orig_dist))

            results = process_contacts_and_calculate_distance(contact_data)

            all_results.append({
                'pdb_file': pdb_file,
                'contacts': results
            })

            cmd.delete("all")
    
    # 将结果保存到JSON文件
    with open(output_json, 'w') as f:
        json.dump(all_results, f, indent=4)
    print(f"Results have been saved to {output_json}")
