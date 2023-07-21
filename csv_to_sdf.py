#!/usr/bin/env python3

import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import os


df = pd.read_csv("/Users/alejandroflores/ado/COMPCHEM/structures/a2a2b_prieto.csv") # Reading the .csv 
df = df.reset_index(drop=True)
print(df.columns) 

'''
                                          Code\tSMILES
0    sy1irp-58\tCCCN1C(NC(C(C(OCC)=O)=C1C)C2=CC=CO2)=O
1    sy1irp-61\tO=C(NC(C(C(OCC)=O)=C1C)C2=CC=CO2)N1...
2    267\tO=C(NC(C(C(OCC)=O)=C1C)C2=CC=CO2)N1CC3=C(...
3    sy1lg-08\tCCCN1C(=O)NC(c2ccco2)C(C(=O)OC(C)C)=C1C
4    sy1lg-09\tCC1=C(C(=O)OC(C)C)C(c2ccco2)NC(=O)N1...
..                                                 ...
221  jl-esp-3b-2\tO=C(OCCC)C1=C(C)N=C2N(C(C=CC=C3)=...
222  sy1rpd-436b\tO=C(OCCCC)C1=C(C)N=C2N(C(C=CC=C3)...
223  jl-esp-3d-2\tO=C(OC)C1=C(C)N=C2N(C(C=CC=C3)=C3...
224  jl-esp-3e-2\tO=C(OCCC)C1=C(C)N=C2N(C(C=CC=C3)=...
225  sy1rpd-437b\tO=C(OCCCC)C1=C(C)N=C2N(C(C=CC=C3)...
'''

# We have one columns. For both Code and SMILES. Code;SMILES. 
# We want to separate Code and SMILES and convert them to 2 new columns
# Apparantly Code and Smiles are separated within the same column by \t so we use str.split function 

df[['Number_Name', 'SMILES']] = df['Code\tSMILES'].str.split('\t', expand=True)
df = df.drop(df.columns[0], axis=1)
print(df)

# Now df should look like this: 

''''
     Number_Name                                             SMILES
0      sy1irp-58             CCCN1C(NC(C(C(OCC)=O)=C1C)C2=CC=CO2)=O
1      sy1irp-61    O=C(NC(C(C(OCC)=O)=C1C)C2=CC=CO2)N1CC3=CC=CC=C3
2            267  O=C(NC(C(C(OCC)=O)=C1C)C2=CC=CO2)N1CC3=C(F)C=C...
3       sy1lg-08            CCCN1C(=O)NC(c2ccco2)C(C(=O)OC(C)C)=C1C
4       sy1lg-09      CC1=C(C(=O)OC(C)C)C(c2ccco2)NC(=O)N1Cc1ccccc1
..           ...                                                ...
221  jl-esp-3b-2  O=C(OCCC)C1=C(C)N=C2N(C(C=CC=C3)=C3N2CC4=C(F)C...
222  sy1rpd-436b  O=C(OCCCC)C1=C(C)N=C2N(C(C=CC=C3)=C3N2CC4=C(F)...
223  jl-esp-3d-2  O=C(OC)C1=C(C)N=C2N(C(C=CC=C3)=C3N2CC4=C(F)C=C...
224  jl-esp-3e-2  O=C(OCCC)C1=C(C)N=C2N(C(C=CC=C3)=C3N2CC4=C(F)C...
225  sy1rpd-437b  O=C(OCCCC)C1=C(C)N=C2N(C(C=CC=C3)=C3N2CC4=C(F)...

[226 rows x 2 columns]
'''

def sdf_from_csv(csv_file):
    #Create a directory to store the individual .sdf files
    output_directory = "Molecules_from_csv"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        
    #Export each molecule to a separate .sdf file
    for i, row in df.iterrows():
        molecule_name = row["Number_Name"]
        smiles_code = row["SMILES"]

        # Convert SMILES to RDKit molecule object
        molecule = Chem.MolFromSmiles(smiles_code)
        if pd.isna(molecule_name) or pd.isna(smiles_code):         # Skip rows with missing values
            print(f"Could not read molecule from row {i + 2}.")
            continue

        # File name 
        filename = os.path.join(output_directory, f"{molecule_name}.sdf")

        # Write molecule to new .sdf file
        writer = Chem.SDWriter(filename)
        writer.write(molecule)
        writer.close()
        
    print(f"All molecules exported to individual .sdf files in the '{output_directory}' directory.")



if __name__ == "__main__":
    # Call the function to convert SMILES to .sdf files
    parser = argparse.ArgumentParser(description="Convert CSV to SDF")
    parser.add_argument("csv_filename", help="The name of the CSV file to convert")
    args = parser.parse_args()
    sdf_from_csv(df)

