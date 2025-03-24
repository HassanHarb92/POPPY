from rdkit import Chem

def atomic_number_to_symbol(atomic_number):
    element = Chem.rdchem.GetPeriodicTable().GetElementSymbol(atomic_number)
    return element

# Example usage:
print(atomic_number_to_symbol(1))  # Should print H
print(atomic_number_to_symbol(6))  # Should print C
print(atomic_number_to_symbol(58))
print(atomic_number_to_symbol(23))
