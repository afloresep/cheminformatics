# Cheminformatics Scripts

This repository contains a collection of cheminformatics scripts written in Python using libraries such as RDKit and pandas. These scripts are designed to assist with various tasks related to chemical structure manipulation, analysis, and data processing.

## List of Scripts

1. `csv_to_sdf.py`: A script to convert a CSV file with molecule names and SMILES codes into separate .sdf files, each named after the molecule.

2. `sdf_reader.py`: A script to convert several molecules in a single .sdf file into separate .sdf file for each molecule


## Getting Started

To use these scripts, you need to have Python and the required libraries installed. You can install the dependencies using pip:
pip install pandas
pip install rdkit

## How to Use 
You can use it as any other python script: 
    python csv_to_sdf.py

You can also make the script executable. In you terminal, navigate to the directory containing your 'csv_to_sdf.py' script and run 
`chmod +x csv_to_sdf.py`
then move the script to a directory that's in the system's PATH variable or add its location to the PATH variable so you can call the script from any location in the terminal. 
Then you can call the Python function from the command line as follows: 
`csv_to_sdf.py your_csv_file.csv`


