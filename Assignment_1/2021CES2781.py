#!/usr/bin/env python
# coding: utf-8

# <h1>Assignment 1

# <h3>Q1: Reversing a list

# In[11]:


#following are sample input and corresponding outputs
Input  = [12, 35, 9, 56, 24]
Output = [24, 35, 9, 56, 12]

Input  = [1, 2, 3]
Output = [3, 2, 1]


# In[12]:


#write your solution here to reverse a list

print(Input[::-1])


# <h3>Q2: Swapping elements in lists

# In[13]:


Input  = [12, 35, 9, 56, 24]
Output = [12, 9, 35, 56, 24]


# In[14]:


#write your solution here

def swap(i,j,S):
    S[i],S[j]=S[j],S[i]
    return S
print(swap(1,2,Input))


# Q3: Reversing a string

# In[15]:


input_str  = 'This is CVL757'
output_str = '757LVC si sihT'


# In[16]:


#write your solution here
print(input_str[::-1])


# <h3>Q4: Reversing the order of words in a string

# In[17]:


input_str  = 'This is CVL757'
output_str = 'CVL757 is This'


# In[18]:


#write your solution here
def reverse_str(a:str,sep=' '):
    a=a.split(sep)
    return " ".join(a[::-1])
print(reverse_str(input_str))


# <h3>Q5: Initialise a 3x3 matrix with random values, use numpy

# In[4]:


import numpy as np
random_matrix = np.random.rand(3,3) #write solution here


# In[6]:


np.random.seed(73)
print(random_matrix)


# <h3>Q6: Swap the columns of random_matrix

# In[7]:


inp = np.array([[0.13806091, 0.82474243, 0.15646752],
                [0.50029824, 0.45875794, 0.53401557],
                [0.95397773, 0.79407795, 0.07442586]])

out = np.array([[0.15646752, 0.82474243, 0.13806091],
                [0.53401557, 0.45875794, 0.50029824],
                [0.07442586, 0.79407795, 0.95397773]])


# In[22]:


#write your solution here to swap 1st and last column
inp[:,[0,2]]=inp[:,[2,0]]
print(inp)


# <h3>Q7: Swap the rows of random_matrix

# In[8]:


inp = np.array([[0.13806091, 0.82474243, 0.15646752],
               [0.50029824, 0.45875794, 0.53401557],
               [0.95397773, 0.79407795, 0.07442586]])

out = np.array([[0.95397773, 0.79407795, 0.07442586],
               [0.50029824, 0.45875794, 0.53401557],
               [0.13806091, 0.82474243, 0.15646752]])


# In[9]:


#write your solution here to swap 1st and 3rd rows
inp[[0,-1],:]=inp[[-1,0],:]
print(inp)


# <h3>Q8: Plot the sine and cosine curves as shown below

# In[11]:


import matplotlib.pyplot as plt


# In[12]:


x = np.arange(0,2*np.pi,0.01)
#write code here to plot and show the legend also
fig,ax=plt.subplots(1,1)
ax.plot(x,np.sin(x))
ax.plot(x,np.cos(x))
ax.set_xlabel("x-------->")
ax.set_ylabel("f(x)-------->")
ax.legend(["Sin(x)","Cos(x)"])
plt.show()


# <h3>Q9: Plot the following truss</h3>
# Length of element 2, 6 and 9 (between nodes 1 and 3, 3 and 5, and 5 and 9) is 5.<br>
# Length of element 3 and 7 is 7m.
#  

# <img src='q9.svg' alt="Truss Q9" width="500" height="600">

# In[13]:


#write solution here
Nodes=[1,2,3,4,5,6]
Nod_cord=[(0,0),(5,7),(5,0),(10,7),(10,0),(15,0)]
Elements=[(1,2),(1,3),(2,3),(2,4),(2,5),(3,5),(4,5),(4,6),(5,6)]
fig,ax=plt.subplots(1,1)
ax.scatter([i[0] for i in Nod_cord],[i[1] for i in Nod_cord],s=100,color='k')
for i in Elements:
    ax.plot([Nod_cord[i[0]-1][0],Nod_cord[i[1]-1][0]],[Nod_cord[i[0]-1][1],Nod_cord[i[1]-1][1]],color='k')
                                                        
ax.set_xlim(-5,20)
ax.set_ylim(-5,10)


# **Q10: Consider the plane truss shown above. Given E = 200GPa and A = 0.005m2,
# determine and horizontal load of 20kN at node 2. Both node 1 and node 6 have pin supports:**
# 
# 1. the global stiffness matrix for the structure.
# 2. the horizontal and vertical displacements at nodes 2, 3, 4, and 5.
# 3. the horizontal and vertical reactions at nodes 1 and 6.
# 4. the stress in each element.

# In[2]:


#write solution here
import numpy as np

class Node:
    """Node
        -Num= node id
        -Pos: Coordinate Position
        -Disp: [ux,uy], =[None,None] if unknown
        -Restrained[i]==1 if restrained 0 if unrestrained
        -Load Vector=[Px,Py]
        -Association[i]=Structure DOF corresponding to ith local dof
        """
    def __init__(self,num,Pos,Disp=[0.0,0.0],Restrain=[0,0],Load=[0,0],Association=[None,None]) -> None:
        self.Num=id
        self.Pos=np.asarray(Pos)
        self.Disp_vec=np.asarray(Disp)
        self.restrained=np.asarray(Restrain)
        self.Load=np.asarray(Load)
        self.Association=Association
    def set_Pos(self,Pos):
        self.Pos=Pos
    def get_Pos(self):
        return self.Pos
    def get_ID(self):

        return self.ID
    def get_Disp(self):
        return self.Disp_vec
    def get_disp_pos(self,exxag=1.0):
        prev_pos=self.Pos.reshape(-1,1)
        Disp_vec=self.Disp_vec.reshape(-1,1)
        new_pos=np.add(prev_pos,exxag*Disp_vec)
        return new_pos
    def set_Disp(self,Disp):
        self.Disp_vec=np.asarray(Disp)
    def get_Load(self):
        return self.Load
    def get_restrains(self):
        return self.restrained
    def set_Load(self,Load=[0,0]):
        self.Load=np.asarray(Load)
    def set_Restr(self,Restrain=[0,0]):
        self.restrained=np.asarray(Restrain)
    def get_Association(self):
        res=[]
        res=self.Association.copy()
        return res
    def set_Association(self,Association :list):
        assert len(Association)==2
        self.Association=Association
class Bar_element:
    """Bar Element
        -E,A,L
        -K[2x2]=[[k -k],[-k,k]]
        -Theta (in degrees) from global axes"""
    def __init__(self,N1:Node,N2 :Node,E=1,A=1) -> None:
        self.L=np.linalg.norm(np.add(N2.get_Pos(),-1*N1.get_Pos()))
        self.N1=N1
        self.N2=N2
        self.E=E
        self.A=A
        self.K=E*A/self.L
        if(np.add(N2.get_Pos(),-1*N1.get_Pos())[0]==0.0):
            self.Theta=90.0
        else:
            self.Theta=np.arctan(np.add(N2.get_Pos(),-1*N1.get_Pos())[1]/np.add(N2.get_Pos(),-1*N1.get_Pos())[0])*180/np.pi
        
    def get_local_k(self):
        #Returns local stiffness matrix [2x2]
        return np.asarray([[self.K,-1*self.K],[-1*self.K,self.K]])

    def get_disp(self):
        return np.concatenate((self.N1.get_Disp(),self.N2.get_Disp()))
    
    def get_force(self):
        return np.concatenate((self.N1.get_Load(),self.N2.get_Load()))
    def get_global_k(self):
        #Return global stiffness matrix =[T]'*[K_Local]*[T]
        #T=[[cos()  sin()   0   0   ],
        #   [0      0   cos()   sin()]]
        if(self.Theta==90.0):
            c=0.0
            s=1.0
        else:
            c=np.cos(self.Theta*np.pi/180)
            s=np.sin(self.Theta*np.pi/180)
        T :np.ndarray=np.asarray([[c,s,0.0,0.0],[0.0,0.0,c,s]])
        K_Local=self.get_local_k()
        return np.dot(np.dot(T.T,K_Local),T)

    def get_stress(self):
        U=self.get_disp()
        #Convert to local coodrinates
        if(self.Theta==90.0):
            c=0.0
            s=1.0
        else:
            c=np.cos(self.Theta*np.pi/180)
            s=np.sin(self.Theta*np.pi/180)
        T :np.ndarray=np.asarray([[c,s,0.0,0.0],[0.0,0.0,c,s]])
        U_local=np.dot(T,U)
        print("Nodal Displacements(Local):",U_local)
        print("Nodal Forces(Local)       :",np.dot(self.get_local_k(),U_local))
        return (U_local[1]-U_local[0])*self.E/self.L
    def set_E(self,E):
        self.E=E
    def set_A(self,A):
        self.A=A
    def get_E(self):
        return self.E
    def get_A(self):
        return self.A
    def get_Nodes(self):
        return [self.N1,self.N2]
    def get_Association(self):
        return self.N1.get_Association()+self.N2.get_Association()

Nod_cord=[(0,0),(5,7),(5,0),(10,7),(10,0),(15,0)]
Elements=[(1,2),(1,3),(2,3),(2,4),(2,5),(3,5),(4,5),(4,6),(5,6)]

Nodes=dict()
Elem=dict()
for i in range(len(Nod_cord)):
    Nodes[i+1]=Node(i+1,[Nod_cord[i][0],Nod_cord[i][1]])
Nodes[2].set_Load([20,0])
Nodes[1].set_Restr([1,1])
Nodes[6].set_Restr([0,1])
for i in range(len(Elements)):
    N1=Elements[i][0]
    N2=Elements[i][1]
    Elem[i+1]=Bar_element(Nodes[N1],Nodes[N2],E=2e8,A=0.005)

#Node numbering with structural dofs, forming association matrices
Count_P=0
DOF_counter=0
for i in Nodes.keys():
    #Take node
    node :Node=Nodes[i]
    restr=node.get_restrains()
    Associate=node.get_Association()
    for j in range(len(restr)):
        if(restr[j]==0):
            Associate[j]=DOF_counter
            Count_P+=1
            DOF_counter+=1
    node.set_Association(Associate)

#Now number restrained joints

for i in Nodes.keys():
    #Take node
    node :Node=Nodes[i]
    #Check its restraints
    restr1=Nodes[i].get_restrains()
    Associate=node.get_Association()
    for j in range(len(restr1)):
        if(restr1[j]==1):
            Associate[j]=DOF_counter
            DOF_counter+=1
    node.set_Association(Associate)

N=2*len(Nodes)
P=Count_P             #No. of Unrestrained joints for KPP :P x P
K_Ts=np.zeros((N,N)) #Total stiffness matrix

for i in Elem.keys():
    member :Bar_element=Elem[i]
    Association=member.get_Association()
    k_g=member.get_global_k()
    K_Ts[np.ix_(Association,Association)]=np.add(K_Ts[np.ix_(Association,Association)],k_g)

F=np.zeros((N,1))     #Nodal Force vector 
U=np.zeros((N,1))     #Nodal Displacement Vector 
for i in Nodes.keys():
    node :Node =Nodes[i]
    Association=node.get_Association()
    node_disp =node.get_Disp().reshape( -1,1)
    node_force=node.get_Load().reshape( -1,1)
    F[Association]=node_force
    U[Association]=node_disp
    
#Partition into Kpp,Kxp,Kpx,Kxx and Up,Ux,Fp,Fx

Kpp =K_Ts[:P,:P]
Kpx =K_Ts[:P,P:]
Kxp =K_Ts[ P:,:P]
Kxx =K_Ts[ P:,P:]

 
Fp=F[:P]
Fx=F[P:]
Up=U[:P]
Ux=U[P:]

"""Solving force-displacement equations
    [Kpp]{Up}+[Kpx]{Ux}={Fp}
    [Kxp]{Up}+[Kxx]{Ux}={Fx}

    {Up}=[Kpp]A(-l){[Fp]-[Kpx][Ux]}     #Unknown forces
    {Fx}=[Kxp]{Up}+[Kxx]{Ux}            #Unknown Reactions
 """
Up=np.dot(np.linalg.inv(Kpp),np.add(Fp,-1*np.dot(Kpx,Ux)))
Fx=np.add(np.dot(Kxp,Up),np.dot(Kxx,Ux)) 
U=np.concatenate((Up,Ux))
F=np.concatenate((Fp,Fx))

#UPDATING NODAL VALUES
for i in Nodes.keys():
    node :Node =Nodes[i]
    Association=node.get_Association()
    node.set_Disp(U[Association])
    node_force=node.set_Load(F[Association])

for i in Elem.keys():
    member :Bar_element =Elem[i]
    print("Element: ",i)
    stress=member.get_stress()
    print("Stress: ",stress)


# <H1>BONUS QUESTION</H1>
# Plot the defomred shape of truss obtained in Q10

# In[22]:


#write solution here
Nod_cord=[(0,0),(5,7),(5,0),(10,7),(10,0),(15,0)]
New_cord=[]
for i in Nodes.keys():
    
    coord=list(Nodes[i].get_disp_pos(2000.0).flat)
    New_cord+=[(coord[0],coord[1])]
Elements=[(1,2),(1,3),(2,3),(2,4),(2,5),(3,5),(4,5),(4,6),(5,6)]

fig,ax=plt.subplots(1,1)
ax.scatter([i[0] for i in New_cord],[i[1] for i in New_cord],s=100,color='k')

for i in Elements:
    ax.plot([New_cord[i[0]-1][0],New_cord[i[1]-1][0]],[New_cord[i[0]-1][1],New_cord[i[1]-1][1]],color='k',linestyle='--')
    ax.plot([Nod_cord[i[0]-1][0],Nod_cord[i[1]-1][0]],[Nod_cord[i[0]-1][1],Nod_cord[i[1]-1][1]],color='b')


# In[ ]:




