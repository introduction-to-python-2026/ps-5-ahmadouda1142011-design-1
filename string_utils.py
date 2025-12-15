def split_before_each_uppercases(formula):
    if formula == "":
        return []

    parts = []
    start = 0

    for i in range(1, len(formula)):
        if formula[i].isupper():
            part = formula[start:i]
            parts.append(part)
            start = i

    last_part = formula[start:]
    parts.append(last_part)

    return parts




def split_at_first_digit(formula):

    for i in range(len(formula)):
        if formula[i].isdigit():
            j = i
            while j < len(formula) and formula[j].isdigit():
                j += 1
            prefix = formula[:i]
            number_text = formula[i:j]
            number = int(number_text)
            return prefix, number
    return formula, 1

def count_atoms_in_molecule(molecular_formula):

    atoms_dict = {}

    for atom in split_before_each_uppercases(molecular_formula):
        atom_name, atom_count = split_at_first_digit(atom)
        atoms_dict[atom_name] = atom_count


    return atoms_dict



def parse_chemical_reaction(reaction_equation):
    """Takes a reaction equation (string) and returns reactants and products as lists.  
    Example: 'H2 + O2 -> H2O' → (['H2', 'O2'], ['H2O'])"""
    reaction_equation = reaction_equation.replace(" ", "")  # Remove spaces for easier parsing
    reactants, products = reaction_equation.split("->")
    return reactants.split("+"), products.split("+")

def count_atoms_in_reaction(molecules_list):
    """Takes a list of molecular formulas and returns a list of atom count dictionaries.  
    Example: ['H2', 'O2'] → [{'H': 2}, {'O': 2}]"""
    molecules_atoms_count = []
    for molecule in molecules_list:
        molecules_atoms_count.append(count_atoms_in_molecule(molecule))
    return molecules_atoms_count



