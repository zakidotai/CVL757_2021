#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Q.1 Solution:
print("Q.1 Solution")
a = [12, 35, 9, 56, 24]
print(a)
a.reverse()
print(a)
b = [1, 2, 3]
print(b)
b.reverse()
print(b)


# In[ ]:


# Q.2Solution
print("Q.2 Solution")
a=[12,35,9,56,24]
print(a)
p1,p2=1,2
a[p1],a[p2]=a[p2],a[p1]
print(a)


# In[ ]:


# Q.3 Solution
print("Q.3 Solution")
a ='This is CVL757'
print(a)
b=reversed(a)
print(''.join(b))


# In[ ]:


# Q.4 Solution:
print("Q.4 Solution")
a= 'This is CVL757'
print(a)
w=str.split(a)
w.reverse()
print(' '.join(w))


# In[ ]:


# Q.5 Solution:
print("Q.5 Solution")
import numpy as np
a=np.random.randn(3,3)
print(a)


# In[ ]:


# Q. 6 Solution:
print("Q.6 Solution")
print(a)
a[:,[0,2]]= a[:,[2,0]]
print(a)


# In[ ]:


#Q.7 solution
print("Q.7 Solution")
print(a)
a[[0,2],:]=a[[2,0],:]
print(a)


# In[ ]:


#Q.8  Solution
print("Q.8 Solution")
import matplotlib.pyplot as plt
import numpy as np
fig, Q=plt.subplots()
x = np.arange(0,2*np.pi,0.01)
y=np.sin(x)
z=np.cos(x)
Q.plot(x,y, label="sinθ")
Q.plot(x,z, label="cosθ")
plt.legend()
plt.show()


# In[ ]:


# Q.9 Solution
print("Q.9 Solution")
import math
import numpy as np
import matplotlib.pyplot as plt
n=6                                      #no. of nodes
e=9                                      # no. of elements
xco=[0,5,5,10,10,15]                     # x-coordinates of node of nodes 1 to 6 respectively
yco=[0,7,0,7,0,0]                        # y-coordinates of node of nodes 1 to 6 respectively
sn=[1,1,3,2,2,3,5,4,5]                   # start node of element of elements 1 to 9 respectively
en=[2,3,2,4,5,5,4,6,6]                   # end node of element of elements 1 to 9 respectively
x=[]
y=[]
for i in range(e):
    s=sn[i]
    e=en[i]
    x+=[xco[s-1],xco[e-1]]
    y+=[yco[s-1],yco[e-1]]
    plt.plot(x,y,color="black")    


# In[ ]:


# Q.10 Solution
print("Q.10 Solution")
import math
from colorama import Fore, Back, Style
import numpy as np
n=6                                      #no. of nodes
e=9                                      # no. of elements
xco=[0,5000,5000,10000,10000,15000]                     # x-coordinates of node of nodes 1 to 6 respectively
yco=[0,7000,0,7000,0,0]                        # y-coordinates of node of nodes 1 to 6 respectively
A=0.005e6                                # Area of Cross section of truss member in mm
E=200e3                                  # Young's Modulus of Elasticity in N/mm or MPa                             
sn=[1,1,3,2,2,3,5,4,5]                   # start node of element of elements 1 to 9 respectively
en=[2,3,2,4,5,5,4,6,6]                   # end node of element of elements 1 to 9 respectively
len=[] 
elcon=[] 
cosofel=[] 
sinofel=[]

for i in range(e):
    x1=float(xco[sn[i]-1])
    x2=float(xco[en[i]-1])
    y1=float(yco[sn[i]-1])
    y2=float(yco[en[i]-1])
    l=math.sqrt((x2-x1)**2+(y2-y1)**2)
    cons=A*E/l
    cos=(x2-x1)/l
    sin=(y2-y1)/l
    len.append(l)
    elcon.append(cons)
    cosofel.append(cos)
    sinofel.append(sin)

    
esm= []                                                 #element stiffness matrix
for i in range(e):
    cc = float(cosofel[i])**2
    ss = float(sinofel[i])**2
    cs = float(cosofel[i])*float(sinofel[i])
    
    mat = elcon[i]*np.array([[cc, cs, -cc, -cs],
                      [cs, ss, -cs, -ss],
                      [-cc, -cs, cc, cs],
                      [-cs, -ss, cs, ss]])

    esm.append(mat)
       
kg =[]                                     # list of Global stiffness matrix (embedded in size of totalstiffness matrix )
for i in range(e):                    
                      
    ass=[2*sn[i]-1,2*sn[i],2*en[i]-1,2*en[i]]                     # association matrix of element i
                                         
    gmat =np.zeros((n*2,n*2))    
    elmat=esm[i]                  
    for j in range(4):                  
        for k in range(4):              
            a=ass[j]-1                     
            b=ass[k]-1                
            gmat[a,b] = elmat[j,k]         # Global stiffness matrix (embedded in size of totalstiffness matrix )
    kg.append(gmat)                        # list of associated global stiffness matrix

kts=np.zeros((n*2,n*2))
for i in range(e):
    kts=kts+kg[i]                          # total stiffness matrix
    
print("global stiffness matrix for the structure is ")
print(kts)
displist=['u1', 'v1', 'u2', 'v2', 'u3', 'v3', 'u4', 'v4', 'u5', 'v5', 'u6', 'v6']
forcelist=['fx1', 'fy1', 'fx2', 'fy2', 'fx3', 'fy3', 'fx4', 'fy4', 'fx5', 'fy5', 'fx6', 'fy6']

dispmat = np.ones((n*2,1))
for i in range(n*2):
    if i==0 :
        dispmat[i,0]=0
    elif i==1 :
        dispmat[i,0]=0
    elif i==10 :
        dispmat[i,0]=0
    elif i==11 :
        dispmat[i,0]=0    
    else:
        dispmat[i,0]=1

forcemat=np.array([[0],[0],[20e3],[0],[0],[0],[0],[0],[0],[0],[0],[0]])
rlist=[]
for i in range(n*2):
    if dispmat[i,0] == 0:
        rlist.append(i)      


rrkts = np.delete(kts, rlist, 0) 
crkts = np.delete(rrkts, rlist, 1) 
kpp= crkts
rforcemat = np.delete(forcemat, rlist, 0) 
rdispmat = np.delete(dispmat, rlist, 0)
dispresult = np.matmul(np.linalg.inv(kpp), rforcemat)
j=0
for i in range(n*2):
    if dispmat[i,0]==1:
        dispmat[i,0]=dispresult[j,0]
        j=j+1
        
print(Fore.RED+"\n\n\n\nThe horizontal displacements at nodes 2 is %f mm" % dispmat[2,0])  
print(Fore.RED+"\n\n\n\nThe Vertical displacements at nodes 2 is %f mm" % dispmat[3,0])
print(Fore.RED+"\n\n\n\nThe horizontal displacements at nodes 3 is %f mm" % dispmat[4,0])
print(Fore.RED+"\n\n\n\nThe Vertical displacements at nodes 3 is %f mm" % dispmat[5,0])
print(Fore.RED+"\n\n\n\nThe horizontal displacements at nodes 4 is %f mm" % dispmat[6,0])
print(Fore.RED+"\n\n\n\nThe Vertical displacements at nodes 4 is %f mm" % dispmat[7,0])
print(Fore.RED+"\n\n\n\nThe horizontal displacements at nodes 5 is %f mm" % dispmat[8,0])
print(Fore.RED+"\n\n\n\nThe Vertical displacements at nodes 5 is %f mm" % dispmat[9,0])

        
forceresult = np.matmul(kts, dispmat)
fkn1=forceresult/1000
print(Fore.BLUE+"\n\n\n\nthe horizontal  reactions at nodes 1 is %f kN" % fkn1[0,0])
print(Fore.BLUE+"\n\n\n\nthe vertical  reactions at nodes 1 is %f kN" % fkn1[1,0])
print(Fore.BLUE+"\n\n\n\nthe horizontal  reactions at nodes 2 is %f kN" % fkn1[10,0])
print(Fore.BLUE+"\n\n\n\nthe vertical  reactions at nodes 2 is %f kN" % fkn1[11,0])
newxco = []
newyco = []
count = 0
for i in range(n):
    k = xco[i]+dispmat[count,0]
    newxco.append(k)
    count = count+1
    l = yco[i]+dispmat[count,0]
    newyco.append(l)
    count = count+1

    newlenofel = []
for i in range(e):
    a, b = sn[i], en[i]
    x1 = float(newxco[a-1])
    y1 = float(newyco[a-1])
    x2 = float(newxco[b-1])
    y2 = float(newyco[b-1])
    l = math.sqrt((x2-x1)**2+(y2-y1)**2)
    newlenofel.append(l)

    
elstrain = np.zeros((e,1))
for i in range(e):
    elstrain[i,0] = (newlenofel[i]-len[i])/(len[i]) 
    
    
elstress = np.zeros((e,1))
for i in range(e):
    elstress[i,0] = E * elstrain[i,0] 
    print(Fore.RED+"\n\nstress in element {0} is {1} Mpa" .format(i+1,elstress[i,0]))
    
    
    


# In[ ]:


# Bonus Question
print("Bonus Question Solution")
import math
from colorama import Fore, Back, Style
import numpy as np
n=6                                      #no. of nodes
e=9                                      # no. of elements
xco=[0,5000,5000,10000,10000,15000]                     # x-coordinates of node of nodes 1 to 6 respectively
yco=[0,7000,0,7000,0,0]                        # y-coordinates of node of nodes 1 to 6 respectively
A=0.0000005e6                                # Area of Cross section of truss member in mm
E=200e3                                  # Young's Modulus of Elasticity in N/mm or MPa                             
sn=[1,1,3,2,2,3,5,4,5]                   # start node of element of elements 1 to 9 respectively
en=[2,3,2,4,5,5,4,6,6]                   # end node of element of elements 1 to 9 respectively
len=[] 
elcon=[] 
cosofel=[] 
sinofel=[]

for i in range(e):
    x1=float(xco[sn[i]-1])
    x2=float(xco[en[i]-1])
    y1=float(yco[sn[i]-1])
    y2=float(yco[en[i]-1])
    l=math.sqrt((x2-x1)**2+(y2-y1)**2)
    cons=A*E/l
    cos=(x2-x1)/l
    sin=(y2-y1)/l
    len.append(l)
    elcon.append(cons)
    cosofel.append(cos)
    sinofel.append(sin)

    
esm= []                                                 #element stiffness matrix
for i in range(e):
    cc = float(cosofel[i])**2
    ss = float(sinofel[i])**2
    cs = float(cosofel[i])*float(sinofel[i])
    
    mat = elcon[i]*np.array([[cc, cs, -cc, -cs],
                      [cs, ss, -cs, -ss],
                      [-cc, -cs, cc, cs],
                      [-cs, -ss, cs, ss]])

    esm.append(mat)
       
kg =[]                                     # list of Global stiffness matrix (embedded in size of totalstiffness matrix )
for i in range(e):                    
                      
    ass=[2*sn[i]-1,2*sn[i],2*en[i]-1,2*en[i]]                     # association matrix of element i
                                         
    gmat =np.zeros((n*2,n*2))    
    elmat=esm[i]                  
    for j in range(4):                  
        for k in range(4):              
            a=ass[j]-1                     
            b=ass[k]-1                
            gmat[a,b] = elmat[j,k]         # Global stiffness matrix (embedded in size of totalstiffness matrix )
    kg.append(gmat)                        # list of associated global stiffness matrix

kts=np.zeros((n*2,n*2))
for i in range(e):
    kts=kts+kg[i]                          # total stiffness matrix
    

displist=['u1', 'v1', 'u2', 'v2', 'u3', 'v3', 'u4', 'v4', 'u5', 'v5', 'u6', 'v6']
forcelist=['fx1', 'fy1', 'fx2', 'fy2', 'fx3', 'fy3', 'fx4', 'fy4', 'fx5', 'fy5', 'fx6', 'fy6']

dispmat = np.ones((n*2,1))
for i in range(n*2):
    if i==0 :
        dispmat[i,0]=0
    elif i==1 :
        dispmat[i,0]=0
    elif i==10 :
        dispmat[i,0]=0
    elif i==11 :
        dispmat[i,0]=0    
    else:
        dispmat[i,0]=1

forcemat=np.array([[0],[0],[20e3],[0],[0],[0],[0],[0],[0],[0],[0],[0]])
rlist=[]
for i in range(n*2):
    if dispmat[i,0] == 0:
        rlist.append(i)      


rrkts = np.delete(kts, rlist, 0) 
crkts = np.delete(rrkts, rlist, 1) 
kpp= crkts
rforcemat = np.delete(forcemat, rlist, 0) 
rdispmat = np.delete(dispmat, rlist, 0)
dispresult = np.matmul(np.linalg.inv(kpp), rforcemat)
j=0
for i in range(n*2):
    if dispmat[i,0]==1:
        dispmat[i,0]=dispresult[j,0]
        j=j+1
        

        
forceresult = np.matmul(kts, dispmat)
fkn1=forceresult/1000

newxco = []
newyco = []
count = 0
for i in range(n):
    k = xco[i]+dispmat[count,0]
    newxco.append(k)
    count = count+1
    l = yco[i]+dispmat[count,0]
    newyco.append(l)
    count = count+1

import math
import numpy as np
import matplotlib.pyplot as plt
n=6                                      #no. of nodes
e=9                                      # no. of elements
xco=[0,5,5,10,10,15]                     # x-coordinates of node of nodes 1 to 6 respectively
yco=[0,7,0,7,0,0]                        # y-coordinates of node of nodes 1 to 6 respectively
sn=[1,1,3,2,2,3,5,4,5]                   # start node of element of elements 1 to 9 respectively
en=[2,3,2,4,5,5,4,6,6]                   # end node of element of elements 1 to 9 respectively
p=[]
q=[]
x=[]
y=[]

for i in range(e):
    s=sn[i]
    e=en[i]
    p+=[newxco[s-1],newxco[e-1]]
    q+=[newyco[s-1],newyco[e-1]]
    plt.plot(p,q,color="darkred")
    
      
   

