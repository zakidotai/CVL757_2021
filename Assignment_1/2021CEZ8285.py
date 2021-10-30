#!/usr/bin/env python
# coding: utf-8

# In[1]:


Input  = [12, 35, 9, 56, 24]
Input.reverse()
print(Input)


# In[2]:


Input  = [12, 35, 9, 56, 24]
v=Input[0]
Input[0]=Input[4]
Input[4]=v
print(Input)


# In[3]:


txt = "This is CVL757"[::-1]
print(txt)


# In[7]:


input = 'This is CVL757'
# first split the string into words 
words = input.split(' ') 
# then reverse the split string list and join using space 
sentence = ' '.join(reversed(words)) 
print(sentence)


# In[19]:


import numpy as np
seed_values=24
np.random.seed(seed_values)
np.random.rand(3,3)


# In[20]:


inp = np.array([[0.13806091, 0.82474243, 0.15646752],
                [0.50029824, 0.45875794, 0.53401557],
                [0.95397773, 0.79407795, 0.07442586]])
inp[:,[2, 0]]=inp[:,[0, 2]]
print(inp)


# In[21]:


inp = np.array([[0.13806091, 0.82474243, 0.15646752],
               [0.50029824, 0.45875794, 0.53401557],
               [0.95397773, 0.79407795, 0.07442586]])
inp[[2, 0],:]=inp[[0, 2],:]
print(inp)


# In[30]:


import matplotlib.pyplot as plt
import numpy as np
x= np.arange(0,2*np.pi,0.01)
y=np.sin(x)
z=np.cos(x)
plt.plot(x,y,x,z)
plt.legend(['sin(\u03B8)', 'cos(\u03B8)'])
plt.show()
plt.savefig('2021CEZ8285(1).png')


# In[31]:


x=[5,5,0,5,10,5,10,10,15,10]
y=[0,7,0,0,0 ,7,7, 0, 0,7]
plt.plot(x,y)
plt.show()
plt.savefig('2021CEZ8285(2).png')


# In[20]:


import numpy as np
A= 0.005
E = 2e8
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
    
print("gstiff")
print(gstiff)

b = np.delete(gstiff,0, 0)
b = np.delete(b,0,1)
b = np.delete(b,0, 0)
b = np.delete(b,0,1)
b = np.delete(b,8, 0)
b = np.delete(b,8,1)
b = np.delete(b,8, 0)
b = np.delete(b,8,1)
f1=np.matrix('20;0;0;0;0;0;0;0')
u = np.linalg.inv(b).dot(f1)
U=np.zeros([tdofs,1])

i=0
while i<8:
    U[i+2][0]=u[i][0]
    i +=1
f= gstiff[0:2].dot(U)
print("reactions at 1st node")
print(f)
f2= gstiff[10:12].dot(U)
print("reactions at 2nd node")
print(f2)
print("nodal displacement")    
print(np.matrix([U[0],U[1],U[10],U[11]]))
print("other displacement")
print(u)
ielem=0

while ielem < 9: 
    x1 =x[ndcon[ielem][0]-1]
    x2 =x[ndcon[ielem][1]-1]
    y1 =y[ndcon[ielem][0]-1]
    y2 =y[ndcon[ielem][1]-1] 
    L1 = mt.sqrt((x2-x1)**2+(y2-y1)**2)
    C =(x2-x1)/L1
    S = (y2-y1)/L1
    stress=(E/L1)*np.matrix([[C,S,-C,-S],[-C,-S,C,S]]).dot(np.matrix([U[ndcon[ielem][0]*2-2],U[ndcon[ielem][0]*2-1],U[ndcon[ielem][1]*2-2],U[ndcon[ielem][1]*2-1]]))
    print("stress in element %s" % (ielem+1)) 
    print(stress)
    ielem +=1


# In[28]:


import numpy as np
import matplotlib.pyplot as plt
ndcon= np.array([[1,2],[1,3],[3, 2],[2, 4],[2, 5],[3, 5],[5, 4],[4, 6],[5,6]])
nelem = len(ndcon)
x=[0,5+2.18759959e-0,5+ 1.11111111e-01,10+ 1.85426626e-0,10+ 2.22222222e-01,15]
y=[0,7+-3.50052910e-01,0-3.50052910e-01,7+ 1.11957672e-01,0+-5.41375661e-01,0]
x1=[0,5,5,10,10,15]
y1=[0,7,0,7,0,0]
def TrussPlot(x,y,line_style):
    ielem = 0
    while ielem < 9: 
        x1 =x[ndcon[ielem][0]-1]
        x2 =x[ndcon[ielem][1]-1]
        y1 =y[ndcon[ielem][0]-1]
        y2 =y[ndcon[ielem][1]-1]
        line_style=str(line_style)
        plt.plot([x1,x2],[y1,y2],line_style)
        ielem+=1
plt.show()

TrussPlot(x,y,'r')
TrussPlot(x1,y1,'b')
plt.savefig('2021CEZ8285(3).png')


# In[ ]:





# In[ ]:




