#coding:utf-8
import numpy as np
from scipy.linalg import sqrtm
from scipy import linalg
import string

f=open('hamsqr1.dat','r')
line=f.readline()
line=f.readline()
line=f.readline()
line=f.readline()
line=f.readline()
line=f.readline()
H=[]
while line:
    line1=line.rstrip('\n')
    a=line1.split(' ')
    c2=[float(s) for s in a if s!='']
    c3=[float(s) for s in c2 if s!='\n']
    n=len(c3)/2
    #=[]
    #H.append(c3)
    for i in range(n):
        hcomp= complex(c3[2*i], c3[(2*i+1)])
        H.append(hcomp)
    #H.append(c2)
    #line=f.readline()
print(H)
mH=np.matrix(H)
#for i in range(8):
#    print(" ")
#    for j in range(8):
#        value = mH[i,j]
#        print("{0} {1} {2}".format(i+1,j+1,(value.real*value.real + value.imag*value.imag)))


#print "H"
#print mH
#print(mH.size)
f.close()

f=open('oversqr.dat', 'r')
line=f.readline()
line=f.readline()
line=f.readline()
line=f.readline()
line=f.readline()
line=f.readline()
S=[]

while line:
    line1=line.rstrip('\n')
    a=line1.split(' ')
    c2=[float(s) for s in a if s!='']
    '''
    c3=[float(s) for s in c2 if s!='\n']
    n=len(c3)/2
    si=[]
    for i in range(n):
        scomp=complex(c3[2*i], c3[2*i+1])
        si.append(scomp)
    '''
    S.append(c2)
    line=f.readline()
#print S
mS=np.matrix(S)
#print "S"
#print mS
#print(mS.size)

f.close()

s2= np.sqrt(2)/2.0
a=0.5
a1=-0.5
#a=complex(s2,0)
#a1=complex(-s2,0)
#b=complex(0,s2)
#b1=complex(0,-s2)
'''
Q=[
   [a,a1,a1,a,0,0,0,0],
   [a,a,a1,a1,0,0,0,0],
   [a,a1,a,a1,0,0,0,0],
   [a,a,a,a,0,0,0,0],
   [0,0,0,0, a,a1,a1,a],
   [0,0,0,0, a,a,a1,a1],
   [0,0,0,0, a,a1,a,a1],
   [0,0,0,0, a,a,a,a]
]
'''

Q=[
   [a,a1,a,a1],
   [a,a1,a1,a],
   [a,a,a1,a1],
   [a,a,a,a]
   ]


#Q=[[1,1,1,1], [1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]]
#Q=[[1,1],[1,-1]]
#Q=[[1,0,0,1,0,0],[0,1,0,0,1,0],[0,0,1,0,0,1], [1,0,0,-1,0,0],[0,1,0,0,-1,0],[0,0,1,0,0,-1]]
#print Q

mQ=np.matrix(Q)
#print("mq")
#print mQ
#mP=np.matrix(P)
#print mP


#print "S^1/2"
#sqrtS=sqrtm(S)
#msqrtS=np.matrix(sqrtS)
#print msqrtS

#print "S^-1/2"
#invsqrtS=linalg.inv(msqrtS)
#print invsqrtS

#print "S^1/2*S^(-1/2)=?I"
#checkinv=msqrtS*np.matrix(invsqrtS)
#print checkinv

#TinvsqrtS=invsqrtS.T
#print "S^-1/2 T"
#print TinvsqrtS

#tmp=np.matrix(TinvsqrtS)*mH
#F=tmp*invsqrtS
#print "S^(-1/2)HS^(1/2)"
#print F

#tmp = F*mQ
#QFQ=np.matrix(mQ.T)*tmp
#print "QFQ"
#print QFQ





#PSP
#SP=mS*mP
#PSP= mP*SP
#print "PSP"
#print PSP


#PHP



#print "S*Q"
SQ=mS*mQ
print(mS)
print(mQ)
print(SQ)

#print SQ

#print "Q^H"
#print mQ.H

#print "QSQ"
QSQ=(mQ.H)*SQ
#print QSQ


#print 'H*Q'
HQ=mH*mQ
#print HQ
#print("Q*H*Q")
QHQ=(mQ.H)*HQ
#print QHQ
for i in range(4):
    print(" ")
    for j in range(4):
        value = QSQ[i,j]
        print("{0} {1} {2}".format(i+1,j+1,(value.real*value.real + value.imag*value.imag)))

#print "(QSQ)^1/2"
r=sqrtm(QSQ)   #r is array , not matrix, we have to do np.matrix(r) transfer it into a matrix, otherwise, the result will be wrong when we do some matrix operations on it
#print r

sqrtQSQ=np.matrix(r)

invsqrtQSQ=linalg.inv(sqrtQSQ)
#print "(QSQ)^(-1/2)"
#print invsqrtQSQ

checkinv=np.matrix(invsqrtQSQ)*np.matrix(r)
#print "check inv"
#print checkinv

tmp=np.matrix(invsqrtQSQ)
HinvsqrtQSQ=tmp.H
FD=QHQ*invsqrtQSQ
DFD=np.matrix(HinvsqrtQSQ)*FD
#print "DFD"
#print DFD

#print "DFD var"
#for i in range(12):
#    print(" ")
#    for j in range(12):
#        value = DFD[i,j]
#        print("{0} {1} {2}".format(i+1,j+1,(value.real*value.real + value.imag*value.imag)))

