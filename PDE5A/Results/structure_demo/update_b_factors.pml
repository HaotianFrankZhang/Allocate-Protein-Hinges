# Load the structure
load 1udu.pdb

# Path to your list of residues and new B-factors
@PDE5A_b_mode1.txt

# Iterate over the file and update B-factors
python
# Open the file
with open('PDE5A_b_mode1.txt', 'r') as file:
    for line in file:
        # Skip comments and empty lines
        if line.startswith('#') or not line.strip():
            continue
        
        # Parse the line
        chain, residue_id, new_b_factor = line.split(',')
        chain = chain.strip()
        residue_id = residue_id.strip()
        new_b_factor = float(new_b_factor.strip())

        # Create a selection for the target residue
        selection = f'/1udu//{chain}/{residue_id}/'
        
        # Update the B-factor
        cmd.alter(selection, f'b={new_b_factor}')
        
        # Rebuild to reflect changes
        cmd.rebuild()
python end

# Save the modified structure
save 1udu_mode1.pdb

