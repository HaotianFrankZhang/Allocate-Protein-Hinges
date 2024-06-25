# Path to your PDB file
pdb_file_path = 'allChemicals_2.pdb'

# Set to store unique ligand names
ligand_names = set()

# Open and read the PDB file
with open(pdb_file_path, 'r') as pdb_file:
    for line in pdb_file:
        # Check if the line corresponds to a heteroatom (non-standard residue, possibly a ligand)
        if line.startswith('HETATM'):
            # Extract the residue name (ligand name) from the line
            # Typically, the ligand name is located between positions 17 and 20 in the PDB file format
            ligand_name = line[17:20].strip()  # Remove any whitespace
            # Add the ligand name to the set (this ensures uniqueness)
            ligand_names.add(ligand_name)

# Print out the unique ligand names
print(ligand_names)

with open('./Moles.txt', 'w') as wf:
    for item in ligand_names:
        wf.write(item + ', ')



