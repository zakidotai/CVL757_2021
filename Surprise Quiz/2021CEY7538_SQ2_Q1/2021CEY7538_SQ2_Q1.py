#!/usr/bin/env python
# coding: utf-8

# ## QUESTION 1
# ### SURPRISE TEST

# #### Material and load values

# In[17]:




E=210e6 #kN/m^2
mu=0.3
t=0.025 #m


# #### STIFFNESS MATRIX OF TRANGULAR ELEMENT

# In[48]:


import numpy as np

def LinearTriangleElementStiffness(E,mu,t,xi,yi,xj,yj,xm,ym,p):
    A = (xi*(yj-ym)+ xj*(ym-yi) + xm*(yi-yj))/2
    bi = yj-ym
    bj = ym-yi
    bm = yi-yj
    gi = xm-xj
    gj = xi-xm
    gm = xj-xi
    B=np.asarray([[bi, 0 ,bj, 0, bm, 0],
                  [0, gi, 0, gj, 0 ,gm],
                  [gi,bi,gj,bj,gm,bm]])

    B=B*(1/(2*A))
    if p==1:    #plane stress condition
        D=(E/(1-mu*mu))*np.asarray([[1,mu,0],[mu,1,0],[0,0,(1-mu)/2]])
    elif p==2: #plane strain condition
        D=(E/((1+mu)*(1-2*mu)))*np.asarray([[1-mu,mu,0],[mu,1-mu,0],[0,0,(1-2*mu)/2]])
    
    return t*A*np.dot(B.T,np.dot(D,B))
    


# In[100]:


#STIFFNEES MATRIXFOR TRAINGLUAR ELEMENT 1
k1=LinearTriangleElementStiffness(E,mu,t,0,0,0.5,0.25,0,0.25,1)
print('--------------K1------------\n')
print(k1)


# In[101]:


#STIFFNEES AMTRIX FOR TRAINGLUAR ELEMENT 1
k2=LinearTriangleElementStiffness(E,mu,t,0,0,0.5,0,0.5,0.25,1)
print('--------------K2------------\n')
print(k2)


# In[71]:


#arrray storing all the traingular element stiffness matrix
km=np.stack((k1,k2),-1)
km=km.transpose(-1,0,1)
print(km.shape)


# ### Global Stiffness Matrix

# In[52]:


#STEP 2 GOLABL MATRIX
K = np.zeros((8, 8))


def LinearTriangleAssemble(K :np.ndarray, k, i, j, m):
    K[2*i-1-1, 2*i-1-1] = K[2*i-1-1, 2*i-1-1] + k[1-1, 1-1]
    K[2*i-1-1, 2*i-1] = K[2*i-1-1,  2*i-1] + k[1-1, 2-1]
    K[2*i-1-1, 2*j-1-1] = K[2*i-1-1,  2*j-1-1] + k[1-1, 3-1]
    K[2*i-1-1,  2*j-1] = K[2*i-1-1,  2*j-1] + k[1-1, 4-1]
    K[2*i-1-1,  2*m-1-1] = K[2*i-1-1,  2*m-1-1] + k[1-1, 5-1]
    K[2*i-1-1,  2*m-1] = K[2*i-1-1,  2*m-1] + k[1-1, 6-1]
    K[2*i-1, 2*i-1-1] = K[2*i-1, 2*i-1-1] + k[2-1, 1-1]
    K[2*i-1, 2*i-1] = K[2*i-1, 2*i-1] + k[2-1, 2-1]
    K[2*i-1, 2*j-1-1] = K[2*i-1, 2*j-1-1] + k[2-1, 3-1]
    K[2*i-1, 2*j-1] = K[2*i-1, 2*j-1] + k[2-1, 4-1]
    K[2*i-1, 2*m-1-1] = K[2*i-1, 2*m-1-1] + k[2-1, 5-1]
    K[2*i-1, 2*m-1] = K[2*i-1, 2*m-1] + k[2-1, 6-1]
    K[2*j-1-1, 2*i-1-1] = K[2*j-1-1,  2*i-1-1] + k[3-1, 1-1]
    K[2*j-1-1,  2*i-1] = K[2*j-1-1,  2*i-1] + k[3-1, 2-1]
    K[2*j-1-1,  2*j-1-1] = K[2*j-1-1,  2*j-1-1] + k[3-1, 3-1]
    K[2*j-1-1,  2*j-1] = K[2*j-1-1,  2*j-1] + k[3-1, 4-1]
    K[2*j-1-1,  2*m-1-1] = K[2*j-1-1,  2*m-1-1] + k[3-1, 5-1]
    K[2*j-1-1,  2*m-1] = K[2*j-1-1, 2*m-1] + k[3-1, 6-1]
    K[2*j-1, 2*i-1-1] = K[2*j-1, 2*i-1-1] + k[4-1, 1-1]
    K[2*j-1, 2*i-1] = K[2*j-1, 2*i-1] + k[4-1, 2-1]
    K[2*j-1, 2*j-1-1] = K[2*j-1, 2*j-1-1] + k[4-1, 3-1]
    K[2*j-1, 2*j-1] = K[2*j-1, 2*j-1] + k[4-1, 4-1]
    K[2*j-1, 2*m-1-1] = K[2*j-1, 2*m-1-1] + k[4-1, 5-1]
    K[2*j-1, 2*m-1] = K[2*j-1, 2*m-1] + k[4-1, 6-1]
    K[2*m-1-1,  2*i-1-1] = K[2*m-1-1,  2*i-1-1] + k[5-1, 1-1]
    K[2*m-1-1,  2*i-1] = K[2*m-1-1, 2*i-1] + k[5-1, 2-1]
    K[2*m-1-1,  2*j-1-1] = K[2*m-1-1,  2*j-1-1] + k[5-1, 3-1]
    K[2*m-1-1,  2*j-1] = K[2*m-1-1,  2*j-1] + k[5-1, 4-1]
    K[2*m-1-1,  2*m-1-1] = K[2*m-1-1,  2*m-1-1] + k[5-1, 5-1]
    K[2*m-1-1, 2*m-1] = K[2*m-1-1,  2*m-1] + k[5-1, 6-1]
    K[2*m-1, 2*i-1-1] = K[2*m-1, 2*i-1-1] + k[6-1, 1-1]
    K[2*m-1, 2*i-1] = K[2*m-1, 2*i-1] + k[6-1, 2-1]
    K[2*m-1, 2*j-1-1] = K[2*m-1, 2*j-1-1] + k[6-1, 3-1]
    K[2*m-1, 2*j-1] = K[2*m-1, 2*j-1] + k[6-1, 4-1]
    K[2*m-1, 2*m-1-1] = K[2*m-1, 2*m-1-1] + k[6-1, 5-1]
    K[2*m-1, 2*m-1] = K[2*m-1, 2*m-1] + k[6-1, 6-1]
    
    return K
    


# In[86]:


K = np.zeros((8, 8))
K=LinearTriangleAssemble(K,km[0],1,3,4)
K=LinearTriangleAssemble(K,km[1],1,2,3)
print("--------GLOBAL MATRIX-------------\n")
print(K)


# ### Boundary conditions

# In[88]:


#Applying Boundary conditions
U=np.zeros((8,1))
F=np.zeros((8,1))

#LOADING CONDITIONS
F[2,0]=9.375
F[3,0]=0
F[4,0]=9.375
F[5,0]=0

#Partitioning THE matrices
Up=U[2:6]
Ux=U[[0,1,6,7],:]
Fp=F[2:6]
Kpp=K[2:6,2:6]
print(Kpp)

print('\n------Known nodal forces-------:=\n')
print(Fp)

Up=np.dot(np.linalg.inv(Kpp),Fp)
print('\n----------Unknown nodal disps----------:=\n')
print(Up)



# ### POST PROCESSING

# In[89]:


#Adding Calcuated unknown displacements to the main displacement vector
U[2:6]=Up

#Calculate force
F=np.dot(K,U)

print('----------FORCES(kN)-------------:=\n')
print(np.round(F,7))


# In[60]:


#Element stresses
def LinearTriangleElementStresses(E,mu,xi,yi,xj,yj,xm,ym,p,u):
    A = (xi*(yj-ym)+ xj*(ym-yi) + xm*(yi-yj))/2
    bi = yj-ym
    bj = ym-yi
    bm = yi-yj
    gi = xm-xj
    gj = xi-xm
    gm = xj-xi
    B=np.asarray([[bi, 0 ,bj, 0, bm, 0],
                  [0, gi, 0, gj, 0 ,gm],
                  [gi,bi,gj,bj,gm,bm]])

    B=B*(1/(2*A))
    
    if p==1:    #plane stress condition
        D=(E/(1-mu*mu))*np.asarray([[1,mu,0],
                                    [mu,1,0],
                                    [0,0,(1-mu)/2]])
        
    elif p==2: #plane strain condition
        D=(E/((1+mu)*(1-2*mu)))*np.asarray([[1-mu,mu,0],
                                            [mu,1-mu,0],
                                            [0,0,(1-2*mu)/2]])
        
    return np.dot(D,np.dot(B,u))


# In[61]:


# FUNCTION TO CALCUATE PRINCIPLE STRESSSES
def LinearTriangleElementPStresses(sigma):
    R = (sigma[0] + sigma[1])/2
    Q = ((sigma[0] - sigma[1])/2)**2 + sigma[2]*sigma[2]
    M = 2*sigma[2]/(sigma[0] - sigma[1])
    s1 = R + np.sqrt(Q)
    s2 = R - np.sqrt(Q)
    theta = (np.arctan(M)/2)*180/np.pi
    return np.asarray([s1 , s2 , theta])


# In[62]:


##ELEMENT STRESSES

#Element disp vector
u1=np.asarray([U[0],U[1],U[4],U[5],U[6],U[7]])
u2=np.asarray([U[0],U[1],U[2],U[3],U[4],U[5]])


# #### Stresses:= [σ<sub>x </sub> σ<sub>y </sub> 𝜏<sub>xy </sub>]

# In[90]:


sigma1=LinearTriangleElementStresses(E,mu,0,0,0.5,0.25,0,0.25,1,u1)
print("-------------Stress(MPa):ELEMENT1-------------\n")
print(np.round(sigma1,6))


# In[91]:


sigma2=LinearTriangleElementStresses(E,mu,0,0,0.5,0,0.5,0.25,1,u2)
print("---------------Stress(MPa):ELEMENT2---------\n")
print(sigma2)


# In[92]:


#Principal stresses
s1=LinearTriangleElementPStresses(sigma1)
print("----------------Principal Stress(MPa):ELEMENT1--------\n")
print(s1)


# In[93]:


s2=LinearTriangleElementPStresses(sigma2)
print("-----------------Principal Stress(MPa):ELEMENT2-------------\n")
print(s2)


# 

# In[ ]:




