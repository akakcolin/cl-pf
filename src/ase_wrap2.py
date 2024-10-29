import os
import pathlib
import argparse
import ase
from ase.io import extxyz
from ase.db import connect
from ase.calculators.singlepoint import SinglePointCalculator


def read_images_from_db(filename):
    db = connect(filename)
    return [image for image in db.select()]

def read_images(filename):
    structs = ase.io.iread(filename)
    return [image for image in structs]

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in', '-i', help='input file type',
                    choices=['vasp', 'cif', 'extxyz', 'asedb'], required=True)
 
    parser.add_argument('--out', '-o', help='output file type',
                    choices=['vasp', 'cif', 'extxyz'], required=True)
    parser.add_argument('input', help='input filename')
    parser.add_argument('output', help='output filename')

    args = parser.parse_args()
    out_filetype = args.out
    in_filetype = args.in
    in_filename = args.input
    out_prex = args.output


    if(in_filetype == 'asedb'):
        images = read_images_from_db(args.input)
        systems_path = args.input.strip('.db')
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
    else:
        images = read_images(args.input)
        index=0
        for a in images:
            outname=out_prex+'_'+str(index)+'.vasp'
            a.write(outname, format=out_filetype, sort=True)
            index=index+1
