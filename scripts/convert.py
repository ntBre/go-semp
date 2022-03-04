#!/usr/bin/python

from sys import argv

PTABLE = {
    "1": "H",
    "6": "C",
}

DERIVED = {
    "DD2",
    "DD3",
    "=",
    "PO1",
    "PO2",
    "PO3",
    "PO7",
    "PO9",
    "EISOL",
    "CORE",
    "EHEAT",
}


def convert(filename):
    """convert the parameters in MOPAC output file `filename` to the format
    expected by MOPAC input files and print the result to stdout

    """
    blank_last = False
    inparam = False
    skip = 0
    with open(filename, "r") as f:
        for line in f:
            sp = line.split()
            if "PARAMETER VALUES" in line:
                skip = 4
                inparam = True
            elif skip > 0:
                skip -= 1
            elif inparam and line == "\n" and blank_last:
                return
            elif inparam and line == "\n":
                blank_last = True
            elif inparam and len(sp) > 3:
                if "XFAC_" in line or "ALPB_" in line:
                    sp[1] = sp[1] + PTABLE[sp[2]]
                    sp[2] = sp[3]
                if sp[1] in DERIVED:
                    continue
                print("%-8s%4s%20.12f" % (sp[1], PTABLE[sp[0]], float(sp[2])))


if __name__ == "__main__":
    if len(argv) < 2:
        print("not enough arguments in call to %s" % argv[0])
        exit(1)
    convert(argv[1])
