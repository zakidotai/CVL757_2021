#!/usr/bin/env python
# coding: utf-8

# # Question : 1

# In[1]:


import numpy as np


# In[2]:


#Given data : material and load 

E=2.1e8 #kN/m^2
v=0.3 
t=0.025 #m
w=100 #kN/m
P=12.5 #kN


# In[ ]:


#Step1- Descritizing the domain
#Step-1 Descritizing the domain
nodes = np.array[(0,0),(0.25,0),(0.125,0.125),(0,0.25),(0.25,0.25),(0.125,0.375),(0,0.5),(0.25,0.5),(0.5,0.25),(0.5,0.5)]
element_nodes = [(1,3,2),(1,4,3),(3,5,2),(3,4,5),(4,6,5),(4,7,6),(5,6,8),(6,7,8),(5,8,9),(8,11,9),(9,11,10)]

def LinearTriangleElementStiffness(E,v,t,xi,yi,xj,yj,xm,ym,p):
    A = (xi*(yj-ym)+ xj*(ym-yi) + xm*(yi-yj))/2
    bi = yj-ym
    bj = ym-yi
    bm = yi-yj
    gi = xm-xj
    gj = xi-xm
    gm = xj-xi
    B=np.asarray([[bi, 0 ,bj, 0, bm, 0],[0, gi, 0, gj, 0 ,gm],[gi,bi,gj,bj,gm,bm]])
    
if p == 1
D = (E/(1-NU*NU))*[1 NU 0 ; NU 1 0 ; 0 0 (1–NU)/2];
elseif p == 2
D = (E/(1+NU)/(1-2*NU))*[1-NU NU 0 ; NU 1–NU 0 ; 0 0 (1–2*NU)/2];
end
y = t*A*B’*D*B

k1=LinearTriangleElementStiffness(E,v,t,0,0,0.5,0.25,0,0.25,1)
print(k1)
k2=LinearTriangleElementStiffness(E,v,t,0,0,0.5,0,0.5,0.25,1)
print(k2)


# In[ ]:




