import os
import pathlib
import argparse
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator

from ase.db import connect
from ase.io import extxyz

def read_images(filename):
    db = connect(filename)
    return [image for image in db.select()]

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('db_file', type=str,
                        help='the .db file to use')

    parser.add_argument("--output", dest="output", type=str, default="out_dftb.xyz",
                        help="output extxyz format")

    args = parser.parse_args()
    outfile = args.output
    systems_path = args.db_file.strip('.db')

    images = read_images(args.db_file)

    systems_dict = {}
    append=True
    for a in images:
        struct = a.toatoms()
        print(struct)
        calculator = SinglePointCalculator(struct,
            energy=a.energy,
            forces=a.forces,
            )
        struct.calc = calculator
        struct.get_potential_energy()
        struct.get_forces()

        extxyz.write_extxyz(outfile, struct, append=append)
