import os
import sys
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt

def ave2sec(x):
  if ( '-' in x ):
    vals = x.split('-')
    times = vals[1].split(':')
    sec = 24*3600*int(vals[0])+3600*int(times[0])+60*int(times[1])+int(times[2])
  else:
    times = x.split(':')
    sec = 3600*int(times[0])+60*int(times[1])+int(times[2])
  return (sec)

f=open("p.dat",'r')
data=[]
for line in f:
    x=line.strip()
    if ('T' in x):
        vals=x.split('T')

        times=vals[1].split(':')
        ymd = vals[0].split('-')
        data.append(vals[0])
timedata=[]
f2=open("p2.dat",'r')
for line in f2:
    x=line.strip()
    timedata.append(ave2sec(x))
#
numcore=[]
f3=open("p3.dat",'r')
for line in f3:
    x=line.strip()
    numcore.append(int(x))
total_time=[]
for i in range(len(numcore)):
    total_time.append(numcore[i]*timedata[i]/3600.0)

print(max(total_time))

#print(i)
a = {}
b={}
index=0
for i in range(len(numcore)):
    d=data[i]
    if d not in a:
        a[d]=1
        b[d]=total_time[i]
    else:
        a[d]+=1
        b[d]+=total_time[i]

total_all=[]
for i in b:
    total_all.append(b[i])
print(max(total_all))
#print(b)
#for i in b:
    #print(b[i])
#a = sorted(a.items(), key=lambda item:item[0])
#print(a)
#result=Counter(data)
total=0.0
for i in b:
    total+=b[i]
print(total/len(b))
labels, values = zip(*Counter(data).items())
labels, values=zip(*b.items())
indexes = np.arange(len(labels))
width = 1
#print(result)
plt.bar(indexes, values, width)
plt.xticks(indexes + width * 1.5, labels)
plt.show()
