from string_utils import split_by_capitals, split_at_number, parse_chemical_reaction, count_atoms_in_molecule, count_atoms_in_reaction
from equation_utils import ELEMENTS, generate_equation_for_element, build_equations, my_solve
from sympy import Eq, symbols, solve

ELEMENTS = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'I', 'Te', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Uut', 'Uuq', 'Uup', 'Uuh', 'Uus', 'Uuo'
]


def generate_equation_for_element(compounds, coefficients, element):
    equation = 0
    for i, compound in enumerate(compounds):
        if element in compound:
            equation += coefficients[i] * compound[element]
    return equation

def build_equations(reactant_atoms, product_atoms):
    reactant_coefficients = list(symbols(f'a0:{len(reactant_atoms)}'))
    product_coefficients = list(symbols(f'b0:{len(product_atoms)}'))
    product_coefficients = product_coefficients[:-1] + [1]

    equations = []
    for element in ELEMENTS:
        lhs = generate_equation_for_element(reactant_atoms, reactant_coefficients, element)
        rhs = generate_equation_for_element(product_atoms, product_coefficients, element)
        if lhs != 0 or rhs != 0:
            equations.append(Eq(lhs, rhs))

    return equations, reactant_coefficients + product_coefficients[:-1]



def my_solve(equations, coefficients):
    solution = solve(equations, coefficients)

    if len(solution) == len(coefficients):
        coefficient_values = []
        for coefficient in coefficients:
            coefficient_values.append(float(solution[coefficient]))
        return coefficient_values


def split_by_capitals(formula):
    t = []
    splited_formula = []
    for char in formula:
        if char.isupper():
            if t:
                splited_formula.append(''.join(t))
            t = [char]
        else:
            t.append(char)
    if t:
        splited_formula.append(''.join(t))
    return splited_formula


def split_at_number(formula):
    digits = []
    letters = []
    for char in formula:
        if char.isdigit():
            digits.append(char)
        else:
            letters.append(char)
    if digits == []:
        digits.append("1")
    return ''.join(letters), int(''.join(digits))

def parse_chemical_reaction(reaction_equation):
    reaction_equation = reaction_equation.replace(" ", "")
    reactants, products = reaction_equation.split("->")
    return reactants.split("+"), products.split("+")


def count_atoms_in_molecule(molecular_formula):
    atom_counts = {}
    for atom in split_by_capitals(molecular_formula):
        atom_name, atom_count = split_at_number(atom)
        atom_counts[atom_name] = atom_counts.get(atom_name, 0) + atom_count
    return atom_counts


def count_atoms_in_reaction(molecules_list):
    molecules_atoms_count = []
    for molecule in molecules_list:
        molecules_atoms_count.append(count_atoms_in_molecule(molecule))
    return molecules_atoms_count
