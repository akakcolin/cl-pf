import os
import sys
import glob
from ccdc.crystal import PackingSimilarity
ie = PackingSimilarity()
from ccdc.io import CrystalReader
c_cellopt2 = [CrystalReader(p)[0] for p in glob.glob("cellopt2_result/*.cif")]                    
c_vasp = [CrystalReader(p)[0] for p in glob.glob("VASP*.cif")]                    
c_cg = [CrystalReader(p)[0] for p in glob.glob("cg_result/*.cif")]                    

s_cellopt2 = [p.split("/")[-1] for p in glob.glob("cellopt2_result/*.cif")]               
s_cg= [p.split("/")[-1] for p in glob.glob("cg_result/*.cif")]               
s_vasp= [p.split("/")[-1] for p in glob.glob("VASP*.cif")]               

for i in range(len(c_cg)) :
    for j in range(len(c_vasp)):
        h=ie.compare(c_cg[i], c_vasp[j])
        print("{0} {1} {2} out of 15 {3}".format(s_cg[i], s_vasp[j], h.nmatched_molecules, h.rmsd))
