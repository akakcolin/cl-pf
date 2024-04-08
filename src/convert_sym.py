from spglib import *
import numpy as np

h1=np.array([[1,0,0], [0,1,0], [0,0,1]])
h2=np.array([[1,0,0], [0,-1,0], [0,0,-1]])
h3=np.array([[-1,0,0], [0,1,0], [0,0,-1]])
h4=np.array([[-1,0,0], [0,-1,0], [0,0,1]])
h5=np.array([[0,1,0], [0,0,1], [1,0,0]])
h6=np.array([[0,1,0], [0,0,-1], [-1,0,0]])
h7=np.array([[0,-1,0], [0,0,1], [-1,0,0]])
h8=np.array([[0,-1,0], [0,0,-1], [1,0,0]])
h9=np.array([[0,0,1], [1,0,0], [0,1,0]])
h10=np.array([[0,0,1], [-1,0,0], [0,-1,0]])
h11=np.array([[0,0,-1], [1,0,0], [0,-1,0]])
h12=np.array([[0,0,-1], [-1,0,0], [0,1,0]])
h13=np.array([[0,-1,0], [-1,0,0], [0,0,-1]])
h14=np.array([[0,-1,0], [1,0,0], [0,0,1]])
h15=np.array([[0,1,0], [-1,0,0], [0,0,1]])
h16=np.array([[0,1,0], [1,0,0], [0,0,-1]])
h17=np.array([[-1,0,0], [0,0,-1], [0,-1,0]])
h18=np.array([[-1,0,0], [0,0,1], [0,1,0]])
h19=np.array([[1,0,0], [0,0,-1], [0,1,0]])
h20=np.array([[1,0,0], [0,0,1], [0,-1,0]])
h21=np.array([[0,0,-1], [0,-1,0], [-1,0,0]])
h22=np.array([[0,0,-1], [0,1,0], [1,0,0]])
h23=np.array([[0,0,-1], [0,-1,0], [1,0,0]])
h24=np.array([[0,0,1], [0,1,0], [-1,0,0]])
h25=np.array([[-1,0,0], [0,-1,0], [0,0,-1]])
h26=np.array([[-1,0,0], [0,1,0], [0,0,1]])
h27=np.array([[1,0,0], [0,-1,0], [0,0,1]])
h28=np.array([[1,0,0], [0,1,0], [0,0,-1]])
h29=np.array([[0,-1,0], [0,0,-1], [-1,0,0]])
h30=np.array([[0,-1,0], [0,0,1], [1,0,0]])
h31=np.array([[0,1,0], [0,0,-1], [1,0,0]])
h32=np.array([[0,1,0], [0,0,1], [-1,0,0]])
h33=np.array([[0,0,-1], [-1,0,0], [0,-1,0]])
h34=np.array([[0,0,-1], [1,0,0], [0,1,0]])
h35=np.array([[0,0,1], [-1,0,0], [0,1,0]])
h36=np.array([[0,0,1], [1,0,0], [0,-1,0]])
h37=np.array([[0,1,0], [1,0,0], [0,0,1]])
h38=np.array([[0,1,0], [-1,0,0], [0,0,-1]])
h39=np.array([[0,-1,0], [1,0,0], [0,0,-1]])
h40=np.array([[0,-1,0], [-1,0,0], [0,0,1]])
h41=np.array([[1,0,0], [0,0,1], [0,1,0]])
h42=np.array([[1,0,0], [0,0,-1], [0,-1,0]])
h43=np.array([[-1,0,0], [0,0,1], [0,-1,0]])
h44=np.array([[-1,0,0], [0,0,-1], [0,1,0]])
h45=np.array([[0,0,1], [0,1,0], [1,0,0]])
h46=np.array([[0,0,1], [0,-1,0], [-1,0,0]])
h47=np.array([[0,0,-1], [0,1,0], [-1,0,0]])
h48=np.array([[0,0,-1], [0,-1,0], [1,0,0]])

#print(symm['rotations'])
a=[i+1 for i in range(48)]
#b=['h'+str(i+1) for i in range(48)]
#print(b)
#mystr = ', '.join(b)
#print( "[ "+ mystr + "]")
h=[ h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, h14, h15, h16, h17, h18, h19, h20, h21, h22, h23, h24, h25, h26, h27, h28, h29, h30, h31, h32, h33, h34, h35, h36, h37, h38, h39, h40, h41, h42, h43, h44, h45, h46, h47, h48]


for i in range(7):
    symm=get_symmetry_from_database(i+1)
    print(symm)
    hsp=get_spacegroup_type(i+1)
    print(hsp)
    hs=get_spacegroup_type(i+1)['hall_symbol']
    trans=symm['translations']
    num=[]
    for r in symm['rotations']:
        for j in zip(h, a):
            if((j[0]==r).all()):
                num.append(j[1])
    #print(num)
    hstr=''
    for n in num:
        hstr=hstr + str(n) + ','

    tstr=''
    for n in trans:
        tstr = tstr + str(n) + ','


    hstr = hstr + " &      !  " + str(i+1) + "    " + hs + "    " + tstr
    print(hstr)

