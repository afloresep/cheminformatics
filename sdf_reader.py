#!bin/bash

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDConfig
import os

def split_sdf_to_files(input_file):
    # Read the .sdf file
    supplier = Chem.SDMolSupplier(input_file)

    # Create a directory to store the individual files
    output_directory = "ligands_exported"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Export each molecule to a separate file
    for i, molecule in enumerate(supplier):
        if molecule is None:
            print(f"Warning: Could not read molecule {i + 1}")
            continue

        # Generate a unique filename for each molecule
        filename = os.path.join(output_directory, f"ligand_{i + 1}.sdf")

        # Write the molecule to a new .sdf file
        writer = Chem.SDWriter(filename)
        writer.write(molecule)
        writer.close()

    print(f"All {i + 1} molecules exported to individual files in the '{output_directory}' directory.")


if __name__ == "__main__":
    input_sdf_file = "/Users/alejandroflores/ado/COMPCHEM/structures/a2a2b_prieto.sdf"  # Replace with the path to your .sdf file
    split_sdf_to_files(input_sdf_file)
