from ase.io import Trajectory, read, write, iread
import argparse

from ase.atoms import *
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import extxyz
import ase.db
import sys
import os
import numpy as np
from ase.units import Bohr, Hartree

if __name__ == "__main__":
    usage = """
    Usage:   geometries.xyz
    stupid code to sort and select xyz structures
    Input Files:
       geometries.xyz   -   xyz-file with geometries, the atom order has to be the same
    Output Files:
       out_sorted.xyz       -   extended xyz-file with geometries, energies and forces for each geometry.

    """
    parser = argparse.ArgumentParser(prog="gather_xyz.py", usage=usage)
    parser.add_argument("geometries_file", type=str, default="geometries.xyz",
                        help="XYZ file containing geometries [default: geometries.xyz]")
    parser.add_argument("sort_file", type=str, default="geometries.xyz",
                        help="XYZ file containing geometries [default: geometries.xyz]")

    parser.add_argument("--output", dest="output", type=str, default="out_sorted.extxyz",
                        help="output extxyz format")

    args = parser.parse_args()

    sort_file = args.sort_file
    outfile = args.output
    sort_energy=[]
    with open(sort_file, 'r') as f:
        line =f.read().split()
        #print(line)
        for l in line:
            sort_energy.append(float(l[7:]))
    print(sort_energy)
    sort_energy = sorted(set(sort_energy))
    sort_energy = list(sort_energy)
    print("length ", len(sort_energy))
    print(sort_energy)

    sort_geo=[]
    for se in sort_energy[:50]:
        structs = ase.io.iread(args.geometries_file)
        for geo in structs:
            #print(geo.info)
            if(abs(geo.info['energy'] - se) < 0.0001):
                print(se)
                sort_geo.append(geo)
                break

    append =True
    # Save results to file as they keep coming in.
    for igeom, results in enumerate(sort_geo):
        extxyz.write_extxyz(outfile, results, append=append)
