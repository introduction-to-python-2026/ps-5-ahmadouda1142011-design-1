from sympy import solve as sympy_solve
from sympy import symbols, Eq
from string_utils import count_atoms_in_molecule

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
    solution = sympy_solve(equations, coefficients)
    if len(solution) == len(coefficients):
        coefficient_values = []
        for c in coefficients:
            coefficient_values.append(solution[c])
        return coefficient_values


def balance_reaction(reaction):
    reactants_str, products_str = reaction.split("->")
    reactants = [r.strip() for r in reactants_str.split("+")]
    products = [p.strip() for p in products_str.split("+")]

    reactant_atoms = [count_atoms_in_molecule(r) for r in reactants]
    product_atoms = [count_atoms_in_molecule(p) for p in products]

    equations, coefficients = build_equations(reactant_atoms, product_atoms)
    values = my_solve(equations, coefficients)

    if values is None:
        return None

    
    return values + [1]
