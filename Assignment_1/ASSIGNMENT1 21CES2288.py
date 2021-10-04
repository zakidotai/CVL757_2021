#!/usr/bin/env python
# coding: utf-8

# <h1>Assignment 1

# <h3>Q1: Reversing a list

# In[9]:


#following are sample input and corresponding outputs
Input  = [12, 35, 9, 56, 24]
Output = [24, 35, 9, 56, 12]

Input  = [1, 2, 3]
Output = [3, 2, 1]


# In[1]:


#write your solution here to reverse a list

input=[1,2,3]
input.reverse()
print(input)


# <h3>Q2: Swapping elements in lists

# In[12]:


Input  = [12, 35, 9, 56, 24]
Output = [12, 9, 35, 56, 24]


# In[2]:


#write your solution here

list = [12, 35, 9, 56, 24]
a = list.index(min(list))
list[a], list[1] = list[1], list[a]
print(list)


# <h3>Q3: Reversing a string

# In[ ]:


input_str  = 'This is CVL757'
output_str = '757LVC si sihT'


# In[3]:


#write your solution here

input_str  = 'This is CVL757'
print(input_str[::-1])


# <h3>Q4: Reversing the order of words in a string

# In[14]:


input_str  = 'This is CVL757'
output_str = 'CVL757 is This'


# In[4]:


#write your solution here

def reverse_str(a:str,sep=' '):
    a=a.split(sep)
    return " ".join(a[::-1])
print(reverse_str(input_str))


# <h3>Q5: Initialise a 3x3 matrix with random values, use numpy

# In[6]:


#write solution here

import numpy as np
x = np.random.random((3,3))
print(x)


# In[23]:





# <h3>Q6: Swap the columns of random_matrix

# In[69]:


inp = np.array([[0.13806091, 0.82474243, 0.15646752],
                [0.50029824, 0.45875794, 0.53401557],
                [0.95397773, 0.79407795, 0.07442586]])

out = np.array([[0.15646752, 0.82474243, 0.13806091],
                [0.53401557, 0.45875794, 0.50029824],
                [0.07442586, 0.79407795, 0.95397773]])


# In[7]:


#write your solution here to swap 1st and last column
a = np.array([[0.365488179, 0.754682,  0.84567182],
              [0.52486712, 0.5468728, 0.54685274],
              [0.78963147, 0.15479635, 0.2846791]])
print(a[:, [2, 1, 0]])


# <h3>Q7: Swap the rows of random_matrix

# In[72]:


inp = np.array([[0.13806091, 0.82474243, 0.15646752],
               [0.50029824, 0.45875794, 0.53401557],
               [0.95397773, 0.79407795, 0.07442586]])

out = np.array([[0.95397773, 0.79407795, 0.07442586],
               [0.50029824, 0.45875794, 0.53401557],
               [0.13806091, 0.82474243, 0.15646752]])


# In[8]:


#write your solution here to swap 1st and 3rd rows
inp = np.array([[0.4586954, 0.84567551,  0.54581552],
              [0.754952142, 0.5452155, 0.45795454],
              [0.55715584, 0.4596327, 0.98524378]])
inp[[0,-1],:]=inp[[-1,0],:]
print(inp)


# <h3>Q8: Plot the sine and cosine curves as shown below

# In[76]:


import matplotlib.pyplot as plt


# In[9]:


x = np.arange(0,2*np.pi,0.01)
#write code here to plot and show the legend also
import matplotlib.pyplot as plt
import numpy as np
x = np.arange(0,2*np.pi,0.01)   # start,stop,step
y = np.sin(x)
z = np.cos(x)
plt.plot(x,y,x,z)
plt.show()


# <h3>Q9: Plot the following truss</h3>
# Length of element 2, 6 and 9 (between nodes 1 and 3, 3 and 5, and 5 and 9) is 5.<br>
# Length of element 3 and 7 is 7m.
#  

# <img src='q9.svg' alt="Truss Q9" width="500" height="600">

# In[10]:


#write solution here
Nodes=[1,2,3,4,5,6]
Nod_cord=[(0,0),(5,7),(5,0),(10,7),(10,0),(15,0)]
Elements=[(1,2),(1,3),(2,3),(2,4),(2,5),(3,5),(4,5),(4,6),(5,6)]
fig,ax=plt.subplots(1,1)
ax.scatter([i[0] for i in Nod_cord],[i[1] for i in Nod_cord],s=100,color='Navy')
for i in Elements:
    ax.plot([Nod_cord[i[0]-1][0],Nod_cord[i[1]-1][0]],[Nod_cord[i[0]-1][1],Nod_cord[i[1]-1][1]],color='Navy')
                                                        
ax.set_xlim(-5,20)
ax.set_ylim(-5,10)


# **Q10: Consider the plane truss shown above. Given E = 200GPa and A = 0.005m2,
# determine and horizontal load of 20kN at node 2. Both node 1 and node 6 have pin supports:**
# 
# 1. the global stiffness matrix for the structure.
# 2. the horizontal and vertical displacements at nodes 2, 3, 4, and 5.
# 3. the horizontal and vertical reactions at nodes 1 and 6.
# 4. the stress in each element.

# In[11]:


#write solution here
#creates a list of member lengths
#the indicies of the list represents the member number as decided

def lmember(coordinates, connections):
    n_cod= coordinates
    conec=connections
    l_mem=[]
    for i in range(len(conec)):
        x1=n_cod [conec[i][0]-1][0]
        x2=n_cod [conec[i][1]-1][0]
        y1=n_cod [conec[i][0]-1][1]
        y2=n_cod [conec[i][1]-1][1]
        length_member=np.sqrt( (x2-x1)**2 +(y2-y1)**2 )
        l_mem.append(length_member)
    return l_mem


#nodes as numbered
coordinates= [(0,0), (5,7), (5,0), (10,7), (10,0), (15,0) ]   #nodes 1,2,3,4,5,6 

#each item represents the node number of the ends of that memeber
connections =[(1,2), (1,3), (3,2), (2,4), (5,2), (3,5), (5,4), (6,4),(5,6)]



c=np.sqrt(5*5 +7*7)
cos=5/c
sin=7/c

#(cos,sin) for each member as per the member connections required for formulating the member stiffness matrix Km
lamxlamy= np.array([ [cos,sin], [1,0], [0,1], [1,0], [-cos,sin], [1,0], [0,1], [-cos,sin], [1,0] ])

#total degrees of freedom = total nodes *2
ndof= len(coordinates)*2
free_dof=8    #these are total known dofs.
dof_seq=[2,3,4,5,6,7,8,9,0,1,10,11]  #write the ID of dofs such that the they are divided into unknown dofs and known dofs
                                    #here, since free_dofs =8, therefore first 8 in dof_seq represent
                                    #known displacements (u-known) and rest are unknown

#Cross sectional properties
E=2e+9
A=0.005


#caclualte the length of members 
L=lmember(coordinates, connections)

#Initialize the member stiffness matricies into a single variable
#Km [member stiiffness] will be of shape ----number of members *4 *4
km=np.zeros(((len(L)),4,4))

for i in range(len(L)):
    c,d =lamxlamy[i]
    km[i]= np.array([[c**2, c*d , -c**2, -c*d],[c*d, d**2, -c*d, -d**2],[-c**2, -c*d, c**2, c*d],[-c*d, -d**2, c*d, d**2 ]])
    gamma= (A*E)/L[i]
    km[i]=km[i] * gamma




#to calcuate the dofs associated with each member.
#dof-label list stores the dofs of each member
#eg. for the first member, the dofs are u0,u1,u2,u2 stored as [0,1,2,3] in the dof_label array


dof_label=[]
dof=np.array([ [1,2], [3,4], [5,6],[7,8], [9,10], [11,12] ])
for p in range(len(connections)):
    i,j =connections[p]
    mm=np.concatenate([dof[i-1],dof[j-1]])-1
    lb=mm.tolist()
    dof_label.append(lb)


##GLOBAL STIFFNESS MATRIX KG
    ##initialized to the size of ndof*ndof

KG= np.zeros((ndof,ndof))

##By using the member stiffness matrix km, the global stiffness matrix is calucated as K
for k in range(km.shape[0]):
    for i in range(4):
        for j in range(4):
            x=int(dof_label[k][j])
            y=int(dof_label[k][i])         
            KG[x,y]= (KG[x,y]+km[k,i,j])
#Printing KG
print(f'Global matrix   \n {KG}')

### MEMBER FORCES AND NODAL DISPLACEMENTS

##Member forces P and nodal displacements as U            
P=np.zeros((ndof,))
U=np.zeros((ndof,))

#Intital condtion of 20KN load along the u3 dof. 

P[2]=20e3

##Shuffle the KG matrix and P such that it is divided into known DOfs and Unknown DOfs in terms of nodal displacements
# the shuffle is done as per dof_seq which specifies the the nodal displacements  "U:= [U_unknows, U_knowns]""
KG[:,:]=KG[dof_seq,:]
KG[:,:]=KG[:,dof_seq]

P[:,]=P[[2,3,4,5,6,7,8,9,0,1,10,11],]

#packing the P as P_knowns(pk) and P_unknowns(pu) [P]:=[pk,pu].  Same for U [U]:=[uu,uk] 
xx=free_dof
pk=P[:xx,]
pu=P[xx:,]
uk=U[xx:,]  
uu=U[:xx,]

#Breaking the shuffled matrix as per uu and uk
K11=KG[:xx,:xx]
K12=KG[:xx,xx:]
K21=KG[xx:,:xx]
K22=KG[xx:,xx:]


#calculate uu(U-unknowns). uu= inverse(K11)*pk
uu= np.matmul(np.linalg.inv(K11),pk)
pu=np.matmul(K21,uu)

#Updating the U and P arrays with calculted values
UU=np.concatenate((uu,uk), axis=0)
FF=np.round(np.concatenate((pk,pu), axis=0),6)
FF2=np.round(np.concatenate((pk,pu), axis=0),6)

#shuffling the U and P to their original sequence as [U]=[U0,U1,U2,U3....] and [P]=[P1,P2,P3......]
UU[2:10,]=UU[:8,]
UU[0]=0
UU[1]=0
node_disp=np.vstack(UU)


FF[2:10,]=FF2[:8,]
FF[0]=FF2[-4]
FF[1]=FF2[-3]
nodal_forces= np.vstack(FF)



#This function calculates the member forces for each member
def memforce(u,m=0, A=0.005,E=2e9):
    u=np.transpose(u)
    lambdas=lamxlamy[m]
    c=lambdas[0]
    d=lambdas[1]
    F=np.zeros((2,))
    t_matrix= np.zeros((2,4))
    alpha=(A*E)/L[m]
    t_matrix= alpha* np.array([ [c,d,-c,-d],[-c,-d,c,d]])
    F= np.matmul(t_matrix,u)
    return F

forces= np.zeros( (2,(len(L))) )
for i in range(len(L)):
    forces[:,i]= np.round(memforce(UU[dof_label[i]] ,m=i),5)
    member_forces= np.vstack(np.array(forces[0]))

###printing member forces, nodal forces and nodal dispalcements
print('\n')
print(f'NODAL DISPLACEMENTS in mm  \n {np.round(node_disp*1000,5)}')
print('\n')
print(f'NODAL FORCES in kN  \n {np.round(nodal_forces/1000,3)}')
print('\n')
print(f'MEMBER FORCES in kN  \n {np.round(member_forces/1000,3)}')
print('\n')
print(f'MEMBER STRESSES in kN/m2  \n {member_forces/(A*1000)}')


# <H1>BONUS QUESTION</H1>
# Plot the defomred shape of truss obtained in Q10

# In[15]:


#write solution here
 
from matplotlib.pyplot import figure
figure(figsize=(8, 6), dpi=60)

##separating nodal displacements into horizontal and vertical displacements
nx=node_disp[::2]
ny=node_disp[1::2]
newcords=[]

#since the displacments are small, they won't be discernable if plotted as per orginal dimensions.fac amplifies the 
#displacements for the sake of plotting
fac=30


##Getting new nodal coordinates bu adjusting displacements
for i in range(len(coordinates)):
    newcords.append((coordinates[i][0]+fac*nx[i][0],coordinates[i][1]+fac*ny[i][0]))

#plotting deformed and undefromed shapes
n_cod=coordinates
conec=connections
n_cod2=newcords
x=[]
y=[]
xx=[]
yy=[]
for i in range(len(conec)):
    x1,x2=n_cod2 [conec[i][0]-1][0],n_cod2 [conec[i][1]-1][0]
    y1,y2=n_cod2 [conec[i][0]-1][1],n_cod2 [conec[i][1]-1][1]
    x.append(x1), x.append(x2),y.append(y1),y.append(y2)
    
    xx1,xx2=n_cod [conec[i][0]-1][0],n_cod [conec[i][1]-1][0]
    yy1,yy2=n_cod [conec[i][0]-1][1],n_cod [conec[i][1]-1][1]
    xx.append(xx1),xx.append(xx2),yy.append(yy1),yy.append(yy2)
plt.plot(xx,yy,'-bo',linewidth=0.5,label='UNDEFORMED')
plt.plot(x,y,'--ro',linewidth=1.0,label='DEFORMED')
plt.legend()
plt.savefig("2021CES2288_QBONUS")
plt.show()


# In[ ]:




