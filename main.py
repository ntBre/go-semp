# dig at levmar
# https://en.wikipedia.org/wiki/Newton%27s_method_in_optimization

# optimize step size in gradient descent
# https://en.wikipedia.org/wiki/Line_search

# gradient descent itself
# https://en.wikipedia.org/wiki/Gradient_descent

from dataclasses import dataclass
import numpy as np
import io

# parameters to skip optimizing because they are either constants or
# derived from other parameters. Gaussian verified to run without
# these
DERIVED_PARAMS = ["PQN", "NValence", "DDN", "KON", "EISol"]

# this says step size should be cube root of machine eps, and that
# machine eps is around 2.2e16, giving 6e-6 for the step
# https://en.wikipedia.org/wiki/Numerical_differentiation
STEP_SIZE = 6e-6

@dataclass
class Atom:
    labl: str
    keys: list[str]
    vals: list[float]

    def values(self) -> np.ndarray:
        return np.array(self.vals, dtype=np.float64)
        
# okay J is definitely a matrix, but there is only one calculation per
# row. f(xi, B) = yi and Ji = df(xi, B)/dB. so each row is the change
# in the function value for a given change in the params. so you do
# finite differences as shown below

# each row has its own value of B, basically vary one component of B
# at a time, central finite differences means

# Ji = df(xi, B)/dB = [ f(xi, B+h) - f(xi, B-h) ] / 2h

# adding in the j index confuses me but I think it will look like this:

# Jij = df(xi, Bj)/dBj = [ f(xi, Bj+h) - f(xi, Bj-h) ] / 2h

# no, that has to be wrong because the beta value to change is fixed
# per row, it doesnt vary with each column

# so there are i rows, and i ranges from 0 to Nparams-1 so there are
# Nparams rows and nstruct columns, which is the transpose of what I
# did last time as expected. so that gives a size of

# n is the number of parameters, so this is backwards and I did it
# right last time

#   m         n

# nparams x nstructs in J and

#    n         m
# nstruct x nparam in Jt so Jt*J is
# nstruct x nstruct and y and f(B) are nstruct x 1

# that doesn"t work! that"s how I got to the transposed form last
# time! Jt is nstruct x nparam and you have to multiply by nstruct x 1
# to get an nparam x

# Jij = what I said above, so each column (j) has a set B or
# equivalently there is one column for each parameter, just like I did
# last time so J is m x n where n is nparams and m is nstruct. then Jt
# is n x m, giving nxn in the product and now y is m x 1 or nstruct x
# 1 as expected

# now the difference between gradient descent and levmar is that
# instead of solving the equation

# (Jt*J)d = Jt*[y - f(B)]

# for d and adding that to B, you just do

# B" = -g * Jt * f(B)

# where g is gamma or the step size, I think

def load_params(filename: str) -> list[Atom]:
    """load semi-empirical parameters from the Gaussian output file
    named by filename and return them as a list of Atoms"""
    atom_label = False
    atom_param = False
    atom = ""
    atoms = []
    a = -1
    with open(filename, "r") as f:
        for line in f:
            if line == " ****\n":
                atom_label = True
            elif atom_label and line == " \n":
                atom_label = False
                atom_param = False
            elif atom_label:
                atom = str(line.split()[0])
                atoms.append(Atom(labl=atom, keys=[], vals=[]))
                atom_label = False
                atom_param = True
                a += 1
            elif atom_param:
                params = line.split()
                for p in params:
                    p, v = p.split("=")
                    if p not in DERIVED_PARAMS:
                        fields = v.split(",")
                        if p == "DCore":
                            p += "=" + ",".join(fields[0:2])
                            fields = fields[2:]
                        for f in fields:
                            atoms[a].keys.append(p)
                            atoms[a].vals.append(f)
    return atoms


def format_params(atoms: list[Atom]) -> str:
    """format semi-empirical parameters of atoms for writing"""
    ret = io.StringIO()
    for atom in atoms:
        ret.write(f"{atom.labl}\n")
        last = None
        for i, k in enumerate(atom.keys):
            if k != last:
                if i > 0:
                    ret.write("\n")
                if "DCore" not in k:
                    ret.write(f"{k}={atom.vals[i]}")
                else:
                    ret.write(f"{k},{atom.vals[i]}")
                last = k
            else:
                ret.write(f",{atom.vals[i]}")
        ret.write("\n****\n")
    return ret.getvalue()


# TODO make this write a real gaussian input file, so it"s going to
# need more args
def dump_params(atoms: list[Atom], filename: str):
    """dump semi-empirical parameters to filename"""
    with open(filename, "w") as out:
        out.write(format_params(atoms))


def load_file07(filename: str):
    """load a list of structures from file07"""
    geoms = []
    with open(filename, "r") as inp:
        buf = []
        for line in inp:
            if "#" not in line:
                fields = line.split()
                buf.extend([float(x) for x in fields])
            elif len(buf) > 0:
                geoms.append(buf)
                buf = []
        geoms.append(buf)
    return geoms


def load_energies(filename: str) -> np.array:
    """load relative ab initio energies from filename"""
    ret = []
    with open(filename, "r") as inp:
        for line in inp:
            ret.append(float(line))
    return np.array(ret, dtype=np.float64)

# TODO write all the gaussian files and corresponding pbs files for a
# set of params


B0 = load_params("opt.out")
dump_params(B0, "test.out")
g = load_file07("file07")
e = load_energies("rel.dat")
