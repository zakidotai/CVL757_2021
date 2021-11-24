#!/usr/bin/env python
# coding: utf-8

# # Surprise TEST 
# ### QUESTION 2 
# 
# 

# #### <font color=black> ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </font>

# #### Material and load values 
# 

# In[3]:


import numpy as np

E=210e6 #kN/m^2
mu=0.3
t=0.025 #m
w=3000 #kN/m^2


# #### STIFFNESS MATRIX OF TRANGULAR ELEMENT 

# In[67]:


## Coordinates information

##coordinates of nodes 
Cords=np.array([[0,0],[0.5,0],[0.5,0.25],[0,0.25],[0.25,0.125]])

#nodes of each triangular element in anticlockwise fashion
Conec=np.array([[1,5,4],[5,3,4],[5,2,3],[1,2,5]])-1     

#number of traingular elements
nEle=Conec.shape[0]

#total number of DOfs
ndof=Cords.shape[0]*2

#dofs corresponding to each node
dofs=np.array([(1,2),(3,4),(5,6),(7,8),(9,10)])-1

#packaging the DOFs as seqence of known DOFs and Unknown Dofs.. dof_seq=[U_unknown U_known]
dof_seq=np.array([1,2,7,8,3,4,5,6,9,10])-1

#total number of free dof's
free_dof=6


# In[68]:


##EleCords contains the coordinates of each node of trangular element
## a typical item represents (x,y) of ith, jth and mth node

EleCords=np.ndarray((nEle, 3,2))  #a typical item represents (x,y) of three nodes of that triangular element in anticlockwise
for i in range(Conec.shape[0]):
    for j in range(Conec.shape[1]):
        EleCords[i,j,0]=Cords[Conec[i][j]][0]
        EleCords[i,j,1]=Cords[Conec[i][j]][1]


# In[69]:


print(EleCords)


# In[74]:


###Cacculates the individual element stiffness matrix of a traingular element

def LinearTriangleElementStiffness(E,mu,t,p,EleNumber=1):
    
    """LinearTraiangleElementStiffnes: E=Elasticity modulus, mu=poissions ratio, t=thickness, 
                                    p=[0 for plane stress, 1 for plane strain]
                                    EleNumber= element number for which the stifness is to be calcuated"""
    
    xi,yi,xj,yj,xm,ym=EleCords[i].flatten()
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
    


# In[75]:


#calcaute individual element stiffness matrix and store it in single array km
km=np.ndarray((nEle,6,6))
for i in range(nEle):
    km[i]=LinearTriangleElementStiffness(E,mu,t,1,EleNumber=i)
print(f'10e4*\n {np.round(km/10000,2)}')


# 
# 

# ### Global Stiffness Matrix 
# 
# ------

# In[76]:


# GOLABL MATRIX
K = np.zeros((8, 8))
#Function calcuates GLOBAL Stiffness matrix

def LinearTriangleAssemble(K :np.ndarray, k, nodes):
    i,j,m=nodes
    
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
    


# In[10]:


np.set_printoptions(linewidth=np.inf)
KG=np.zeros((ndof,ndof))
## Calcuating the GLOBAL MATRIX using individual matricies
for i in range(nEle):
    nodes=Conec[i]+1
    KG=LinearTriangleAssemble(KG,km[i],nodes)
print("--------GLOBAL MATRIX-------------\n")
print(f'10e4*\n {np.round(KG/10000,4)}')   


# -----

# ### Boundary Conditions

# In[77]:


F=np.zeros((ndof,))
U=np.zeros((ndof,))


# In[78]:


#Applying Boundary conditions
U=np.zeros((ndof,1))
F=np.zeros((ndof,1))

#LOADING CONDITIONS
F[2,0]=9.375
F[3,0]=0
F[4,0]=9.375
F[5,0]=0
F[8,0]=0
F[9,0]=0


dof_seq=np.array([3,4,5,6,9,10,1,2,7,8])-1  #segretated as per knows and unknows [knows unknows]
xx=free_dof
#Partitioning THE matrices
Up=U[dof_seq[:xx],]
Ux=U[dof_seq[xx:],] #knows disps
Fp=F[dof_seq[:xx],]
Kpp=KG[np.ix_(dof_seq[:xx],dof_seq[:xx])]


# print(Kpp)

print('\n------Known nodal forces-------:=\n')
print(Fp)

Up=np.dot(np.linalg.inv(Kpp),Fp)
print('\n----------Unknown nodal disps----------:=\n')
print(Up)



# In[81]:


#Adding Calcuated unknown displacements to the main displacement vector
U[dof_seq[:xx],]=Up

#Calculate force
F=np.dot(KG,U)

print('----------FORCES VECTOR(kN)-------------:=\n')
print(np.round(F,7))


# ### POST PROCESSING 

# In[83]:


##Array which contains displacement(nodal) of a traingular element as rows. total colums equals total memebrs


u=np.ndarray((nEle,6))
for p in range(nEle):
    xi=dofs[Conec[p][0]][0]
    yi=dofs[Conec[p][0]][1]
    
    xj=dofs[Conec[p][1]][0]
    yj=dofs[Conec[p][1]][1]
    
    xm=dofs[Conec[p][2]][0]
    ym=dofs[Conec[p][2]][1]
    
    u[p]=U[[xi,yi,xj,yj,xm,ym],].T
    


# In[84]:


# Each row corresponds to displacements of a node of a particular trinagular element
#since each trinagular element has 3 nodes, therefore each row has 6 displacemnts 
#total rows corresponds to total number of memebrs
print(u)


# In[53]:


#Element stresses

def LinearTriangleElementStresses(E,mu,p=1,EleNumber=1):
    xi,yi,xj,yj,xm,ym=EleCords[i].flatten()
    uu=u[EleNumber]
    
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

    return np.dot(D,np.dot(B,uu)).reshape(np.dot(D,np.dot(B,uu)).shape[0],1)


# In[54]:


# FUNCTION TO CALCUATE PRINCIPLE STRESSSES
def LinearTriangleElementPStresses(sigma):
    R = (sigma[0] + sigma[1])/2
    Q = ((sigma[0] - sigma[1])/2)**2 + sigma[2]*sigma[2]
    M = 2*sigma[2]/(sigma[0] - sigma[1])
    s1 = R + np.sqrt(Q)
    s2 = R - np.sqrt(Q)
    theta = (np.arctan(M)/2)*180/np.pi
    return np.asarray([s1 , s2 , theta])


# In[64]:



#sigma is an array having stesses as its items
#for a typical item, it correspongs to sigma_x,sigma_y, tau_xy
#number of items correponds to number of traingular elements

sigma=np.ndarray((nEle,3,1))
for i in range(nEle):
    sigma[i]=LinearTriangleElementStresses(E,mu,EleNumber=i)
# print(sigma*1000)    


# #### Stresses:= [σ<sub>x </sub> σ<sub>y </sub> 𝜏<sub>xy </sub>] 

# In[85]:


for i in range(nEle):
    print(f'\nStress:ELEMENT{i+1}\n\n {sigma[i]}')


# In[86]:


#prinStress is an array having principal stresses and theta items

prinStress=np.ndarray(sigma.shape)
for i in range(nEle):
    prinStress[i]= LinearTriangleElementPStresses(sigma[i])
# print(prinStress)


# ####  Principal Stresses:= [σ<sub>x </sub> σ<sub>y </sub> θ ] 

# In[88]:


for i in range(nEle):
    print(f'\nPrincipal Stress:=ELEMENT{i+1}\n\n {prinStress[i]}')


# In[ ]:





# In[ ]:




