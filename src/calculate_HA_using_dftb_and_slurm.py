#!/usr/bin/env python

import os
import subprocess
import tempfile
import argparse

#
def read_parse_cifs_create_ph_dir(filename):
    input = open(filename, 'r').read().split('\n')
    for i in range(0, len(input)):
        cif_line = input[i]
        dir_line= cif_line.split('_')[0]
        os.mkdir(dir_line)
        os.system('cp ' +str(cif_line) + ' ' + dir_line)
        os.system('cp run_dftb_opt.py  ' + dir_line)
        print("cd " + str(dir_line) + " && python run_dftb_opt.py " + cif_line + "  && cd ..")

# identify opt state
def tail_dynlog(filename):
    input = open(filename, 'r').read().split('\n')
    for i in range(0, len(input)):
        cif_line = input[i]
        #print(line)
        dir_line= cif_line.split('_')[0]
        print("echo "+dir_line)
        #os.systemy("tail -1 "+str(dir_line)+'/dyn.log')
        print("tail -1 "+str(dir_line)+'/dyn.log')
        os.system('cp run_phonopy_dftb.py  ' + dir_line)
        os.system('cp dftb_in.hsd_gamma  ' + dir_line)

def prepare_phonopy_step_one(filename):
    input = open(filename, 'r').read().split('\n')
    for i in range(0, len(input)):
        cif_line = input[i]
        dir_line= cif_line.split('_')[0]
        os.system('cp run_phonopy_dftb.py  ' + dir_line)
        os.system('cp dftb_in.hsd_gamma  ' + dir_line)
        print("cd " + str(dir_line) + " && python run_phonopy_dftb.py geo_end.gen && cd ..")

def run_split_file(splitLen=50, namebase='F', filename='fix_ph_jobs'):
    outputBase =namebase
    input = open(filename, 'r').read().split('\n')
    at = 1
    outlist=[]
    for lines in range(0, len(input), splitLen):
        outputData = input[lines:lines+splitLen]
        newfile=outputBase + str(at).zfill(5)
        outlist.append(newfile)
        output = open(newfile, 'w')
        output.write('\n'.join(outputData))
        output.close()
        at += 1
    for i in outlist:
        dftb_command="parallel -j 4  -a " + str(i)
        gene="python sbatch_submit.py -i \" " + dftb_command + "\" -j " + str(i) +" -f " + str(i)
        os.system(gene)
    for i in outlist:
         print("sbatch run_"+str(i)+".sh")

def prepare_phonopy_step_two(filename):
    input = open(filename, 'r').read().split('\n')
    for i in range(0, len(input)):
        cif_line = input[i]
        dir_line= cif_line.split('_')[0]
        print(" cat " + str(dir_line)+"/run_all_jobs >> all_jobs")

def get_force_phonopy(filename):
    input = open(filename, 'r').read().split('\n')
    for i in range(0, len(input)):
        cif_line = input[i]
        dir_line= cif_line.split('_')[0]
        print("cd " + str(dir_line) + " && phonopy -f disp-???/results.tag --dftb+")

def get_thermal_phonopy(filename):
    input = open(filename, 'r').read().split('\n')
    for i in range(0, len(input)):
        cif_line = input[i]
        dir_line= cif_line.split('_')[0]
        print("cd " + str(dir_line) + " && phonopy -t mesh.conf --dftb+")

def get_zpe_phonopy(filename):
    input = open(filename, 'r').read().split('\n')
    for i in range(0, len(input)):
        cif_line = input[i]
        dir_line= cif_line.split('_')[0]
        print("  grep zero_point_energy "+str(dir_line)+"/thermal_properties.yaml  | awk '{print $2}'")


def get_phonopy_disp_blankdir(filename):
    with open(filename, 'r') as f:
        for line in f:
            li = line.strip()
            lin = str(li).split("&")[0]
            #print(lin)
            dirname = lin.split(' ')[1]
            #print(dirname)
            dftbout = dirname+"/dftb.out"
            if not os.path.isfile(dftbout):
                #print(dirname)
                print(li)
                #newcommand=lin lin + " && cat dftb.out | tail -1"
            #print(newcommand)

def get_phonopy_disp_errdir(filename):
    with open(filename, 'r') as f:
        for line in f:
            li = line.strip()
            lin = str(li).split("&")[0]
            #print(lin)
            dirname = lin.split(' ')[1]
            #print(dirname)
            dftbout = dirname+"/dftb.out"
            if os.path.isfile(dftbout):
                filesize = os.path.getsize(dftbout)
                if filesize == 0:
                    print(li)
                else:
                    with open(dftbout, 'r') as fp:
                        lines = fp.readlines()
                        last_line = lines[-1]
                    if(last_line[0] != '-'):
                        print(li)

def process_phonon_result():
    for istep in range(-5, 9):
        voldir='volume_'+str(istep)
        os.chdir(voldir)
        os.system('phonopy -f disp-???/results.tag --dftb+')
        os.system('phonopy -t mesh.conf --dftb+')
        os.system('cp thermal_properties.yaml ../thermal_properties.yaml_'+str(istep))
        os.chdir('../')

if __name__ == "__main__":
   get_phonopy_disp_errdir('big00001')
   get_zpe_phonopy("cif_file")

