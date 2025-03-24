from rdkit import Chem

def atomic_number_to_symbol(atomic_number):
    element = Chem.rdchem.GetPeriodicTable().GetElementSymbol(atomic_number)
    return element

def extract_initial_structure(log_file_path, xyz_file_path):
    with open(log_file_path, 'r') as file:
        lines = file.readlines()
    
    start_index = 0
    for i, line in enumerate(lines):
        if "Point Number" in line:
            start_index = i
            break

    input_orientation_found = False
    for i in range(start_index, len(lines)):
        if "Input orientation:" in lines[i]:
            input_orientation_found = True
            start_index = i
            break
    
    if not input_orientation_found:
        print("Input orientation not found.")
        return

    # Skip header lines to reach the atomic coordinates
    start_index += 5  # Adjust based on the exact number of header lines before the atomic data

    atoms = []
    while "-----" not in lines[start_index]:
        parts = lines[start_index].split()
        atomic_number = int(parts[1])
        x, y, z = parts[3:6]
        symbol = atomic_number_to_symbol(atomic_number)
        atoms.append((symbol, x, y, z))
        start_index += 1

    with open(xyz_file_path, 'w') as xyz_file:
        xyz_file.write(f"{len(atoms)}\n")
        xyz_file.write("Initial structure\n")
        for atom in atoms:
            xyz_file.write(f"{atom[0]} {atom[1]} {atom[2]} {atom[3]}\n")

# Example usage
log_file_path = 'H-N-IRC.log'  # Update this with your actual log file path
xyz_file_path = 'initial_structure.xyz'
extract_initial_structure(log_file_path, xyz_file_path)

