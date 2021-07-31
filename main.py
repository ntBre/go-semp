from dataclasses import dataclass
import numpy as np
import io
import os

# parameters to skip optimizing because they are either constants or
# derived from other parameters
DERIVED_PARAMS = ['PQN', 'NValence', 'DDN', 'KON', 'EISol']
CHARGE = 0
SPIN = 1


@dataclass
class Atom:
    labl: str
    keys: list[str]
    vals: list[float]

    def values(self) -> np.ndarray:
        return np.array(self.vals, dtype=np.float64)


# just keep duplicating the key !! ex:
# [... Beta Beta Beta Beta ...]
# [... -0.3  0.0  0.0  0.0 ...]

# then when I write them, join with comma while the key is the same

# maybe Atom dataclass with Label and Value fields, then return value
# is a list of Atoms like it is now

# dig at levmar
# https://en.wikipedia.org/wiki/Newton%27s_method_in_optimization

# optimize step size in gradient descent
# https://en.wikipedia.org/wiki/Line_search

# gradient descent itself
# https://en.wikipedia.org/wiki/Gradient_descent

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

# that doesn't work! that's how I got to the transposed form last
# time! Jt is nstruct x nparam and you have to multiply by nstruct x 1
# to get an nparam x

# Jij = what I said above, so each column (j) has a set B or
# equivalently there is one column for each parameter, just like I did
# last time so J is m x n where n is nparams and m is nstruct. then Jt
# is n x m, giving nxn in the product and now y is m x 1 or nstruct x
# 1 as expected

# see this for ideas about step size in qffs and here
# https://en.wikipedia.org/wiki/Numerical_differentiation

# based on that, h should be the cube root of machine epsilon, which
# the article puts around 2.2e-16 for doubles = 6e-6

# now the difference between gradient descent and levmar is that
# instead of solving the equation

# (Jt*J)d = Jt*[y - f(B)]

# for d and adding that to B, you just do

# B' = -g * Jt * f(B)

# where g is gamma or the step size, I think


def load_params(filename: str) -> list[Atom]:
    """load semi-empirical parameters from the Gaussian output file
    named by filename and return them as a list of Atoms"""
    atom_label = False
    atom_param = False
    atom = ""
    atoms = []
    a = -1
    with open(filename, 'r') as f:
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
                    p, v = p.split('=')
                    if p not in DERIVED_PARAMS:
                        fields = v.split(',')
                        if p == 'DCore':
                            p += '=' + ','.join(fields[0:2])
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
                if 'DCore' not in k:
                    ret.write(f"{k}={atom.vals[i]}")
                else:
                    ret.write(f"{k},{atom.vals[i]}")
                last = k
            else:
                ret.write(f",{atom.vals[i]}")
        ret.write("\n****\n")
    return ret.getvalue()


def dump_params(atoms: list[Atom], filename: str):
    """dump semi-empirical parameters to filename in the format expected
    by Gaussian

    """
    with open(filename, 'w') as out:
        out.write(format_params(atoms))


def load_file07(filename: str):
    """load a list of structures from file07"""
    geoms = []
    with open(filename, 'r') as inp:
        buf = []
        for line in inp:
            if '#' not in line:
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
    with open(filename, 'r') as inp:
        for line in inp:
            ret.append(float(line))
    return np.array(ret, dtype=np.float64)


def zip_geom(atoms: list[str], cords: list[float]) -> str:
    """combine a list of atom names with a list of coordinates"""
    ret = io.StringIO()
    for c in range(0, len(cords), 3):
        ret.write("%s%20.12f%20.12f%20.12f\n" %
                  (atoms[c // 3], *cords[c:c+3]))
    return ret.getvalue()


def make_com(comfile: str, paramfile: str,
             atoms: list[str], cords: list[float]):
    """make a Gaussian input file with name comfile, paramfile loaded
    using the @filename construct, and a geometry made by calling
    zip_geom on atoms and cords

    """
    geom = zip_geom(atoms, cords)
    contents = f"""%mem=1000mb
%nprocs=1
#P PM6=(print,zero)

the title

{CHARGE} {SPIN}
{geom}

@{paramfile}

"""
    with open(comfile, 'w') as com:
        com.write(contents)


def make_pbs(pbsfile: str, comfile: str):
    outfile = os.path.splitext(comfile)[0]+'.out'
    contents = f"""#!/bin/sh
#PBS -N sempirical
#PBS -S /bin/bash
#PBS -j oe
#PBS -m abe
#PBS -l mem=1gb
#PBS -l nodes=1:ppn=1

scrdir=/tmp/$USER.$PBS_JOBID

mkdir -p $scrdir
export GAUSS_SCRDIR=$scrdir
export OMP_NUM_THREADS=1

echo "exec_host = $HOSTNAME"

if [[ $HOSTNAME =~ cn([0-9]{{3}}) ]];
then
  nodenum=${{BASH_REMATCH[1]}};
  nodenum=$((10#$nodenum));
  echo $nodenum

  if (( $nodenum <= 29 ))
  then
    echo "Using AVX version";
    export g16root=/usr/local/apps/gaussian/g16-c01-avx/
  elif (( $nodenum > 29 ))
  then
    echo "Using AVX2 version";
    export g16root=/usr/local/apps/gaussian/g16-c01-avx2/
  else
    echo "Unexpected condition!"
    exit 1;
  fi
else
  echo "Not on a compute node!"
  exit 1;
fi

cd $PBS_O_WORKDIR
. $g16root/g16/bsd/g16.profile
g16 {comfile} {outfile}

rm -r $scrdir

"""
    with open(pbsfile, 'w') as pbs:
        pbs.write(contents)


def main():
    labels = ["H", "O", "H"]
    min_geom = 34
    B0 = load_params("opt.out")
    dump_params(B0, "test.out")
    g = load_file07('file07')
    # e = load_energies('rel.dat')
    make_com("test.com", "params.dat", labels, g[min_geom])
    make_pbs("test.pbs", "test.com")


main()
