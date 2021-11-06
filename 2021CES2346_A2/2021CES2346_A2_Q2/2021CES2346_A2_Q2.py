#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math
import numpy as np
te=3
tn=te+1
lem=[3000,4000,2000]
xco=[0,lem[0],(lem[0]+lem[1]),(lem[0]+lem[1]+lem[2])]
yco=[0,0,0,0]

I =5e6
E =210e3

snofel= [1,2,3]                                                   #start node of elements
enofel= [2,3,4]                                                   #end node of elements

lenofel = []                                                  #length of the element
cosofel = []                                                  #cos of element
sinofel = []                                                  #sin of element
elstmat = []                                                  #element stiffness matrix

for i in range(te):  
    a = snofel[i]
    b = enofel[i]
    x1 = float(xco[a-1])
    y1 = float(yco[a-1])
    x2 = float(xco[b-1])
    y2 = float(yco[b-1])
    L= math.sqrt((x2-x1)**2+(y2-y1)**2)
    mat =E*I*np.array([[12/L**3,  6/L**2,  -12/L**3,  6/L**2],
                      [6/L**2  ,     4/L,   -6/L**2,     2/L],
                      [-12/L**3, -6/L**2,   12/L**3, -6/L**2],
                      [6/L**2  ,     2/L,   -6/L**2,    4/L]])
    snofel.append(a)
    enofel.append(b)
    lenofel.append(L)
    elstmat.append(mat)
    
egsmm=[]                               #element global stiffness matrix mapping
for i in range(te):                     
    m = snofel[i]*2                     
    n = enofel[i]*2                     
    assmat = [m-1, m, n-1, n]             # association matrix for element
    
    gmat = np.zeros((tn*2, tn*2))    ## global stiffness matrix loaded with zeros for element(i)
    elmat = elstmat[i]                  ## taking the element stiffness matrix of element(i)
    for j in range(4):                  
        for k in range(4):              
            a = assmat[j]-1                ## addressing row of GST matrix for element(i)
            b = assmat[k]-1                ## addressing column of GST matrix for element(i)
            gmat[a,b] = elmat[j,k]      ## updating the values in GST matrix with EST matrix of element(i)
    egsmm.append(gmat)

GSM = np.zeros((tn*2, tn*2))         ## creating an empyty GSM matrix
for mat in egsmm:
    GSM = GSM+mat  
print("Global stiffness matrix for the structure:\n\n")    
print(GSM)     


displist = []
forcelist = []
for i in range(tn):
    a = str('v')+str(i+1)               # vertical displacement
    displist.append(a)                      
    b = str('w')+str(i+1)               # rotation at node
    displist.append(b) 
    c = str('fy')+str(i+1)              # vertical reaction 
    forcelist.append(c)
    d = str('M')+str(i+1)               # moment at node
    forcelist.append(d) 
dispmat = np.ones((tn*2,1))    

dispmat = np.ones((tn*2,1))
dispmat[1*2-2, 0] = 0
dispmat[2*2-2, 0] = 0
dispmat[3*2-2, 0] = 0
dispmat[4*2-1, 0] = 0
dispmat[4*2-2, 0] = 0

forcemat=np.zeros((tn*2,1))
elnbr=2
W=7
forcemat[(2*elnbr-1),0]=forcemat[(2*elnbr-1),0]-(W*lenofel[elnbr-1]**2)/12
forcemat[(2*elnbr-2),0]= forcemat[(2*elnbr-2),0]-(W*lenofel[elnbr-1])/2  
forcemat[(2*elnbr),0]=forcemat[(2*elnbr),0]-(W*lenofel[elnbr-1])/2 
forcemat[(2*elnbr+1),0]=forcemat[(2*elnbr+1),0]+(W*lenofel[elnbr-1]**2)/12
fcorr=-forcemat


# In[2]:


import numpy
rcdlist = []
for i in range(tn*2):
    if dispmat[i,0] == 0:
        rcdlist.append(i)
rrgsm = numpy.delete(GSM, rcdlist, 0) #row reduction
crgsm = numpy.delete(rrgsm, rcdlist, 1) #column reduction
rgsm = crgsm #reduced global stiffness matrix
rforcemat = numpy.delete(forcemat, rcdlist, 0) #reduced force mat
rdispmat = numpy.delete(dispmat, rcdlist, 0) #reduced disp mat

###_______________Solving____________________###

dispresult = numpy.matmul(numpy.linalg.inv(rgsm), rforcemat)
rin = 0
for i in range(tn*2):
    if dispmat[i,0] == 1:
        dispmat[i,0] = dispresult[rin,0]
        rin = rin+1
#print(dispmat)
forceresult = numpy.matmul(GSM, dispmat)
forceresult=forceresult+fcorr
#print(forceresult)


# In[3]:



print("\n\n\n\nRotations at nodes 1 is  %f rad" % dispmat[3,0])
print("(negative sign indicates clockwise rotation)")
print("\n\n\n\nRotations at nodes 2 is  %f rad" % dispmat[3,0])
print("(negative sign indicates clockwise rotation)")
print("\n\n\n\nRotations at nodes 3 is  %f rad" % dispmat[5,0])


print("\n\n\n\nVertical Reactions at nodes 1 is  %f KN" % (forceresult[0,0]/1000))
print("\n\n\n\nVertical Reactions at nodes 2 is  %f KN" % (forceresult[2,0]/1000))
print("\n\n\n\nVertical Reactions at nodes 3 is  %f KN" % (forceresult[4,0]/1000))
print("\n\n\n\nVertical Reactions at nodes 4 is  %f KN" % (forceresult[6,0]/1000))
print("\n\n\n\nMoment at nodes 4 is  %f KNm" % (forceresult[7,0]/1000000))



u1=np.array([[dispmat[0,0]],[dispmat[1,0]],[dispmat[2,0]],[dispmat[3,0]]])
#print(u1)
u2=np.array([[dispmat[2,0]],[dispmat[3,0]],[dispmat[4,0]],[dispmat[5,0]]])
#print(u2)
u3=np.array([[dispmat[4,0]],[dispmat[5,0]],[dispmat[6,0]],[dispmat[7,0]]])
print("\n\n\n")
f1=np.matmul(elstmat[0],u1)
#print(f1)
f2=np.matmul(elstmat[1],u2)
#print(f2)
f3=np.matmul(elstmat[2],u3)
#print(f3)


# In[4]:


##SHEAR FORCE DIAGRAM
xx=[]
SF=[]
for i in range(1801):
    if i==0:
        V=0
    if i>0 and i/200<3:
        V=forceresult[0]/1000
    if i/200==3:
        V=0
    if i/200>3 and i/200<7:
        V=(forceresult[0]+forceresult[2])/1000-7*(i/200-3)
    if i/200==7:
        V=0    
    if i/200>7 and i/200<9:
        V= -forceresult[6] /1000
    if i/200==9:
        V=0
    xx.append(i/200)    
    SF.append(V) 
  
    
from matplotlib import pyplot as plt
plt.figure(figsize=(16,4))
plt.plot(xx,SF)
plt.plot([0]*10, color='k')
plt.title('SHEAR FORCE DIAGRAM')
plt.xlabel('x m')
plt.ylabel('Shear Force in KN')
plt.ylim(-18,18)
plt.xlim(-1,10)


# In[5]:


## BENDING MOMENT DIAGRAM
xx=[]
BM=[]
for i in range(1801):
    if i==0:
        M=0
    if i>0 and i/200<=3:
        M=forceresult[0]*0.001*i/200
    if i/200>3 and i/200<=7:
        M=(forceresult[0]*i/200+forceresult[2]*(i/200-3))/1000-7*(((i/200-3)**2)/2) 
    if i/200>7 and i/200<9:
        M= forceresult[6]*(9-i/200)/1000+forceresult[7]/1000000
    if i/200==9:
        M=0
    xx.append(i/200)    
    BM.append(M) 
  
    
from matplotlib import pyplot as plt
plt.figure(figsize=(16,4))
plt.plot(xx,BM)
plt.plot([0]*10, color='k')
plt.title('BENDING MOMENT DIAGRAM')
plt.xlabel('x m')
plt.ylabel('BENDING MOMENT in KNm')
#plt.ylim(-18,18)
plt.xlim(-1,10)


# In[ ]:




