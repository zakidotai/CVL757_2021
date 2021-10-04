#!/usr/bin/env python
# coding: utf-8

# In[1]:


# SWAPPING ELEMENTS IN A LIST

def swaplist(num):
 num[0],num[4]=num[4],num[0]

 return num
num=[12,35,9,56,24]
print(swaplist(num))


# In[2]:


# REVERSING A LIST

lst = [1, 2, 3]
lst.reverse()
print(lst)


# In[3]:


# SWAPPING ELEMENTS IN A LIST

def swaplist(num):
 num[1],num[2]=num[2],num[1]

 return num
num=[12,35,9,56,24]
print(swaplist(num))


# In[4]:


# REVERSING A STRING

text = "This is CVL757" [::-1]
print(text)


# In[5]:


# REVERSING ORDER OF WORDS IN A STRING

s  = 'This is CVL757'
l=s.split()
print(l)
l1=l[::-1]
print(l1)
output=' '.join(l1)
print(output)


# In[6]:


# Initialise a 3x3 matrix with random values

import numpy as np

array = np.random.rand(3,3)

print(array)


# In[7]:


# swapping columns of above matrix

array[:, [2,0]] = array[:, [0,2]]
print(array)


# In[8]:


# swapping rows of above matrix

array[[2,0] ,:] = array[[0,2] ,:]
print(array)


# In[14]:


# Plotting the sine and cosine curves

import matplotlib.pyplot as plt
import numpy as np
import math


X = np.arange(0, math.pi*2, 0.05)

y = np.sin(X)
z = np.cos(X)

plt.plot(X, y, color='b', label='sin')
plt.plot(X, z, color='m', label='cos')

plt.xlabel("Angle")
plt.ylabel("Magnitude")
plt.title("Sine and Cosine functions")

plt.legend()

plt.show()


# In[54]:


# plotting the truss

import numpy as np
import matplotlib.pyplot as plt
ndcon= np.array([[1,2],[1,3],[3, 2],[2, 4],[2, 5],[3, 5],[5, 4],[4, 6],[5,6]])
nelem = len(ndcon)
x=[0,5,5,10,10,15]
y=[0,7,0,7,0,0]
ielem = 0
while ielem < 9: 
    x1 =x[ndcon[ielem][0]-1]
    x2 =x[ndcon[ielem][1]-1]
    y1 =y[ndcon[ielem][0]-1]
    y2 =y[ndcon[ielem][1]-1]
    plt.plot([x1,x2],[y1,y2],'Black')
    ielem+=1
plt.show()


# In[1]:


import numpy as np
A= 0.005  # meter^2
E = 2e8   # KN/m^2
L1=5      # meter
L2=7      # meter
L3=8.6023 # meter
ndcon= np.array([[1,2],[1,3],[3, 2],[2, 4],[2, 5],[3, 5],[5, 4],[4, 6],[5,6]])
nelem = len(ndcon)
x=[0,5,5,10,10,15]
y=[0,7,0,7,0,0]
nodes = len(x)
ndofn=2
nnode = 2
tdofs = nodes*ndofn
gstiff = np.zeros([tdofs,tdofs])
import math as mt
ielem=0
while ielem < 9:   
    x1 =x[ndcon[ielem][0]-1]
    x2 =x[ndcon[ielem][1]-1]
    y1 =y[ndcon[ielem][0]-1]
    y2 =y[ndcon[ielem][1]-1]
    gbdof=[]
    L1 = mt.sqrt((x2-x1)**2+(y2-y1)**2)
    C =(x2-x1)/L1
    S = (y2-y1)/L1
    A1 = A*E/L1
    K1 =A1* np.matrix('%s, %s , %s, %s; %s, %s, %s, %s; %s, %s, %s, %s; %s, %s, %s, %s' % (C**2, C*S,-C**2, -C*S, C*S, S**2, -C*S, -S**2,-C**2, -C*S, C**2, C*S, -C*S, -S**2, C*S, S**2))
    
    inode =0
    while inode < nnode:
        idofn=1
        while idofn <= ndofn:
            gbdof.append((ndcon[ielem,inode]-1)*ndofn+idofn)
            idofn += 1
        inode += 1
    i = 0
    if i == 0:
        gstiff[gbdof[i]-1,gbdof[i]-1]=gstiff[gbdof[i]-1,gbdof[i]-1] + K1[0,0]
        gstiff[gbdof[i]-1,gbdof[i+1]-1]=gstiff[gbdof[i]-1,gbdof[i+1]-1] + K1[0,1]
        gstiff[gbdof[i]-1,gbdof[i+2]-1]=gstiff[gbdof[i]-1,gbdof[i+2]-1] + K1[0,2]
        gstiff[gbdof[i]-1,gbdof[i+3]-1]=gstiff[gbdof[i]-1,gbdof[i+3]-1] + K1[0,3]
    i = 1
    if i == 1:
        gstiff[gbdof[i]-1,gbdof[i-1]-1]=gstiff[gbdof[i]-1,gbdof[i-1]-1] + K1[1,0]
        gstiff[gbdof[i]-1,gbdof[i]-1]=gstiff[gbdof[i]-1,gbdof[i]-1]     + K1[1,1]
        gstiff[gbdof[i]-1,gbdof[i+1]-1]=gstiff[gbdof[i]-1,gbdof[i+1]-1] + K1[1,2]
        gstiff[gbdof[i]-1,gbdof[i+2]-1]=gstiff[gbdof[i]-1,gbdof[i+2]-1] + K1[1,3]
    i = 2
    if i == 2:
        gstiff[gbdof[i]-1,gbdof[i-2]-1]=gstiff[gbdof[i]-1,gbdof[i-2]-1] + K1[2,0]
        gstiff[gbdof[i]-1,gbdof[i-1]-1]=gstiff[gbdof[i]-1,gbdof[i-1]-1] + K1[2,1]
        gstiff[gbdof[i]-1,gbdof[i]-1]=gstiff[gbdof[i]-1,gbdof[i]-1]     + K1[2,2]
        gstiff[gbdof[i]-1,gbdof[i+1]-1]=gstiff[gbdof[i]-1,gbdof[i+1]-1] + K1[2,3]
    i = 3 
    if i == 3:
        gstiff[gbdof[i]-1,gbdof[i-3]-1]=gstiff[gbdof[i]-1,gbdof[i-3]-1] + K1[3,0]
        gstiff[gbdof[i]-1,gbdof[i-2]-1]=gstiff[gbdof[i]-1,gbdof[i-2]-1] + K1[3,1]
        gstiff[gbdof[i]-1,gbdof[i-1]-1]=gstiff[gbdof[i]-1,gbdof[i-1]-1] + K1[3,2]
        gstiff[gbdof[i]-1,gbdof[i]-1]=gstiff[gbdof[i]-1,gbdof[i]-1]     + K1[3,3]
    ielem += 1
    
print("Global Stiffness Matrix = ")
print(gstiff)


# In[4]:


b = np.delete(gstiff,0, 0)
b = np.delete(b,0,1)
b = np.delete(b,0, 0)
b = np.delete(b,0,1)
b = np.delete(b,8, 0)
b = np.delete(b,8,1)
b = np.delete(b,8, 0)
b = np.delete(b,8,1)
force=np.matrix('0;0;20;0;0;0;0;0;0;0;0;0')
f1 = np.delete(force,0, 0)
f1 = np.delete(f1,0,0)
f1= np.delete(f1,8, 0)
f1= np.delete(f1,8,0)
u = np.linalg.inv(b).dot(f1)
U=np.zeros([tdofs,1])

i=0
while i<8:
    U[i+2][0]=u[i][0]
    i +=1
print("Displacement at Joints")    
print(np.matrix([U[0],U[1],U[10],U[11]]))
print("other displacement")
print(u)
ielem=0
print("Positive sign = Tensile Force \nNegative Sign = Compression")


# In[5]:


while ielem < 9: 
    x1 =x[ndcon[ielem][0]-1]
    x2 =x[ndcon[ielem][1]-1]
    y1 =y[ndcon[ielem][0]-1]
    y2 =y[ndcon[ielem][1]-1] 
    L1 = mt.sqrt((x2-x1)**2+(y2-y1)**2)
    C =(x2-x1)/L1
    S = (y2-y1)/L1
    stress=(E/L1)*np.matrix([-C,-S,C,S]).dot(np.matrix([U[ndcon[ielem][0]*2-2],U[ndcon[ielem][0]*2-1],U[ndcon[ielem][1]*2-2],U[ndcon[ielem][1]*2-1]]))
    print(f"stress in element {ielem+1} is : {stress} KN/m^2")
    ielem +=1
print("\n Calculation of reactions at node 1 and 6")


# In[7]:


# Calculation of the reactions    
F = np.matrix(gstiff).dot(np.matrix(U))

print(f"Horizontal Reaction at node 1 = {F[0]}\nVertical Reactiona at node 1 = {F[1]}\nHorizontal Reaction at node 6 = {F[10]}\nVertical Reactiona at node 6 = {F[11]}")


# In[ ]:




