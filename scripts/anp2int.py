#!/usr/bin/python3

import re
import numpy as np
from sys import argv

# convert the displacements in an anpass input file to the format expected by
# intder


def parse_anpass(filename):
    """print displacements in filename to intder format and return the
    corresponding energies"""
    # line like "(12F12.8,f20.12)" with optional spaces at each end
    start = re.compile(r"^\s*\(\d+[fF]\d+\.\d,[fF]\d+\.\d+\)\s*$")
    in_disp = False
    energies = []
    count = 0
    with open(filename, "r") as inf:
        for line in inf:
            if start.match(line):
                in_disp = True
            elif "UNKNOWNS" in line:
                print("DISP%5d" % count)
                return energies
            elif in_disp:
                line = [float(x) for x in line.split()]
                energies.append(line.pop())
                for i, disp in enumerate(line):
                    if disp != 0.0:
                        print("%5d%18.8f" % (i, disp))
                print("    0")
                count += 1


def make_rel(energies):
    min_ = min(energies)
    return [e - min_ for e in energies]


if __name__ == "__main__":
    if len(argv) < 2:
        print("not enough arguments in call to %s" % argv[0])
        print("Usage: %s <anpass input file>" % argv[0])
        exit(1)
    energies = make_rel(parse_anpass(argv[1]))
    with open("rel.dat", "w") as out:
        for e in energies:
            out.write("%20.12f\n" % e)
