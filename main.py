# parameters to skip optimizing because they are either constants or
# derived from other parameters
DERIVED_PARAMS = ['PQN', 'NValence', 'DDN', 'KON', 'EISol']

# safe to skip PQN and Nvalence
# same for DDN through EISol

# I really need the params to be ordered to make them easier to put
# into matrices, so two parallel lists with params and the
# corresponding value in the same positions would be ideal, but this
# complicates parsing since some of the values have multiple values
def load_params(filename):
    atom_label = False
    atom_param = False
    atom = ""
    atoms = {}
    with open(filename, 'r') as f:
        for line in f:
            if "parameters" in line:
                print(line, end='')
            elif line == " ****\n":
                atom_label = True
            elif atom_label and line == " \n":
                atom_label = False
                atom_param = False
            elif atom_label:
                atom = str(line.split()[0])
                atoms[atom] = {}
                atom_label = False
                atom_param = True
            elif atom_param:
                params = line.split()
                for p in params:
                    p, v = p.split('=')
                    if p not in DERIVED_PARAMS:
                        atoms[atom][p] = v
    return atoms


print(load_params("opt.out"))
