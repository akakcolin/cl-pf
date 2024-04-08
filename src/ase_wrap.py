#!/usr/bin/env python3
import argparse
parser = argparse.ArgumentParser(description='ase wrapper  ')

parser.add_argument('--inp', '-i', help='input file type',
                    choices=['vasp', 'cif'], required=True)
parser.add_argument('--out', '-o', help='output file type',
                    choices=['vasp', 'cif'], required=True)
parser.add_argument('input', help='input filename')
parser.add_argument('output', help='output filename')


if __name__ == '__main__':
    args = parser.parse_args()
    in_filetype = args.inp
    out_filetype = args.out
    in_filename = args.input
    out_filename = args.output

    # ====================
    #
    print(" Input File                     : %s" % in_filename)
    print(" Output File                    : %s" % out_filename)
    print("")
    #
    from ase import io
    atoms = io.read(in_filename)
    atoms.write(out_filename, format=out_filetype)
