# -*- coding: utf-8 -*-
from ligand_workflow._00_align_ligand import load_and_align_structures
from ligand_workflow._10_calculate_contacts import process_directory as process_contacts
from ligand_workflow._20_convert_json_to_csv import convert_json_to_csv

# ligand_workflow/ligand_workflow.py

def ligand_workflow(input_folder, output_csv, output_json, config_file):
    """
    Main workflow function: align ligands, calculate contacts, and convert JSON to CSV.
    """
    from ._00_super_rnas import super_rnas
    from ._01_ligand_rmsd import ligand_rmsd
    from ._02_find_massC import find_massC
    from ._10_calculate_contacts import process_directory as process_contacts
    from ._20_convert_json_to_csv import convert_json_to_csv

    # Step 1: Align ligands
    print("Aligning ligands...")
    load_and_align_structures(input_folder, config_file)

    # Step 2: Calculate contacts
    print("Calculating contacts...")
    process_contacts(input_folder, config_file, output_json)

    # Step 3: Convert JSON to CSV
    print("Converting JSON to CSV...")
    convert_json_to_csv(output_json, config_file, output_csv)

    print("Workflow completed!")
