import os
import sys

def get_blankdir(filename):
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

def get_errdir(filename):
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
get_errdir('big00001') # file list
