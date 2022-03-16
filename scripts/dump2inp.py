#!/usr/bin/python

import io
from sys import argv

# convert an anpass-style dump file from `pbqff` to rel.dat and file07 expected
# by `semp`

START_GEOM = "# GEOMUP #################"


def format_geom(geom):
    """geom: str -> str"""
    ret = io.StringIO()
    ret.write(START_GEOM + "\n")
    for i in range(0, len(geom), 3):
        for j in range(i, i + 3):
            ret.write("%20.10f" % float(geom[j]))
        ret.write("\n")
    return ret.getvalue()


def work(filename):
    energies = []
    with open(filename, "r") as inp, open("file07", "w") as file07:
        for i, line in enumerate(inp):
            sp = line.split()
            file07.write(format_geom(sp[:-1]))
            energies.append(float(sp[-1]))

    m = min(energies)
    energies = [e - m for e in energies]

    with open("rel.dat", "w") as rel:
        for e in energies:
            rel.write("%20.12f\n" % e)


if __name__ == "__main__":
    if len(argv) < 2:
        print("not enough arguments in call to %s" % argv[0])
        exit(1)
    work(argv[1])
