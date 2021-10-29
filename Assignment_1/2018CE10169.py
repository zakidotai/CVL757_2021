# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% [markdown]
# <h1>Assignment 1
# %% [markdown]
# <h3>Q1: Reversing a list

# %%
#following are sample input and corresponding outputs
Input  = [12, 35, 9, 56, 24]
Output = [24, 35, 9, 56, 12]

#Input  = [1, 2, 3]
#Output = [3, 2, 1]


# %%
#write your solution here to reverse a list

print(Input[::-1])

# %% [markdown]
# <h3>Q2: Swapping elements in lists

# %%
Input  = [12, 35, 9, 56, 24]
Output = [12, 9, 35, 56, 24]


# %%

Input[1],Input[2]=Input[2],Input[1]
print(Input)

# %% [markdown]
# <h3>Q3: Reversing a string

# %%
input_str  = 'This is CVL757'
output_str = '757LVC si sihT'


# %%
#write your solution here
print(input_str[::-1])

# %% [markdown]
# <h3>Q4: Reversing the order of words in a string

# %%
input_str  = 'This is CVL757'
output_str = 'CVL757 is This'


# %%
#write your solution here
def rev_words(a :str,sep=' '):
    a=a.split(sep)
    return " ".join(a[::-1])

print(rev_words(input_str))

# %% [markdown]
# <h3>Q5: Initialise a 3x3 matrix with random values, use numpy

# %%
import numpy as np
        #write solution here


# %%

np.random.seed(50)
random_matrix =np.random.rand(3,3) 
print(random_matrix)

# %% [markdown]
# <h3>Q6: Swap the columns of random_matrix

# %%
inp = np.array([[0.13806091, 0.82474243, 0.15646752],
                [0.50029824, 0.45875794, 0.53401557],
                [0.95397773, 0.79407795, 0.07442586]])

out = np.array([[0.15646752, 0.82474243, 0.13806091],
                [0.53401557, 0.45875794, 0.50029824],
                [0.07442586, 0.79407795, 0.95397773]])

print(inp)


# %%
#write your solution here to swap 1st and last column
inp[:,[0,-1]]=inp[:,[-1,0]]

print(inp)

# %% [markdown]
# <h3>Q7: Swap the rows of random_matrix

# %%
inp = np.array([[0.13806091, 0.82474243, 0.15646752],
               [0.50029824, 0.45875794, 0.53401557],
               [0.95397773, 0.79407795, 0.07442586]])

out = np.array([[0.95397773, 0.79407795, 0.07442586],
               [0.50029824, 0.45875794, 0.53401557],
               [0.13806091, 0.82474243, 0.15646752]])


# %%
#write your solution here to swap 1st and 3rd rows
inp[[0,2],:]=inp[[2,0],:]
print(inp)

# %% [markdown]
# <h3>Q8: Plot the sine and cosine curves as shown below

# %%
import matplotlib.pyplot as plt


# %%
x = np.arange(0,2*np.pi,0.01)
#write code here to plot and_and_and show the legend also
fig,ax=plt.subplots(1,1)
ax.plot(x,np.sin(x))
ax.plot(x,np.cos(x))
ax.set_xlabel("x------->")
ax.set_ylabel("f(x)------>")
ax.legend(["Sin(x)","Cos(x)"])
fig.savefig("2018CE10169_Q8.png")
plt.show()

# %% [markdown]
# <h3>Q9: Plot the following truss</h3>
# Length of element 2, 6 and 9 (between nodes 1 and 3, 3 and 5, and 5 and 9) is 5.<br>
# Length of element 3 and 7 is 7m.
#  
# %% [markdown]
# <img src='q9.svg' alt="Truss Q9" width="500" height="600">

# %%
#write solution here
Nodes=[1,2,3,4,5,6]
Nod_cord=[(0,0),(5,7),(5,0),(10,7),(10,0),(15,0)]
Elements=[(1,2),(1,3),(2,3),(2,4),(2,5),(3,5),(4,5),(4,6),(5,6)]

fig,ax=plt.subplots(1,1)
ax.scatter([i[0] for i in Nod_cord],[i[1] for i in Nod_cord],s=100,color='k')

for i in Elements:
    ax.plot([Nod_cord[i[0]-1][0],Nod_cord[i[1]-1][0]],[Nod_cord[i[0]-1][1],Nod_cord[i[1]-1][1]],color='k')

ax.set_xlim(-2,17)
ax.set_ylim(-2,10)
fig.savefig("2018CE10169_Q9.png")
plt.show()

# %% [markdown]
# **Q10: Consider the plane truss shown above. Given E = 200GPa and A = 0.005m2, and horizontal load of 20kN at node 2. Both node 1 and node 6 have pin supports:**
# 
# 1. the global stiffness matrix for the structure.
# 2. the horizontal and vertical displacements at nodes 2, 3, 4, and 5.
# 3. the horizontal and vertical reactions at nodes 1 and 6.
# 4. the stress in each element.

# %%
#Defining a Node and Bar element
import numpy as np
class Node:
    """Node
        -Num= node id
        -Pos: [x,y] Coordinate Position
        -Disp: [ux,uy],( =[None,None] if unknown)
        -Restrained[i]==1 if restrained 0 if unrestrained
        -Load Vector=[Px,Py]    (Default [0,0])
        -Association :list ->Association[i]=Structure DOF corresponding to ith local DOF
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
    def get_new_pos(self,scale=1.0):
        prev_pos=self.Pos.reshape(-1,1)
        Disp_vec=self.Disp_vec.reshape(-1,1)
        new_pos=np.add(prev_pos,scale*Disp_vec)
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
        -E,A,L  (Material and Geometrical properties)
        -K[2x2]=[[k -k],[-k,k]] (Local Stiffness matrix where k=EA/L)
        -Theta (in degrees) from global x-axis
        -N1, N2 :Nodes of element"""
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
        #print("Nodal Displacements(Local):",U_local)
        #print("Nodal Forces(Local)       :",np.dot(self.get_local_k(),U_local))
        return (U_local[1]-U_local[0])*self.E/self.L
    def get_Nodal_Force(self):
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
        #print("Nodal Displacements(Local):",U_local)
        #print("Nodal Forces(Local)       :",np.dot(self.get_local_k(),U_local))
        return np.dot(self.get_local_k(),U_local)
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
    


# %%
#write solution here 

#units= kN-m

#Step1: Define the Node coordinates and elements as set of two nodes
Nod_cord=[(0,0),(5,7),(5,0),(10,7),(10,0),(15,0)]
Elements=[(1,2),(1,3),(2,3),(2,4),(2,5),(3,5),(4,5),(4,6),(5,6)]

#Step2: Create Nodes and elements
Nodes=dict()
Elem=dict()
for i in range(len(Nod_cord)):
    Nodes[i+1]=Node(i+1,[Nod_cord[i][0],Nod_cord[i][1]])

#Step3: Apply Loads and Restraints
Nodes[2].set_Load([20,0])
Nodes[1].set_Restr([1,1])   #1 for restrained, 0 for unrestrained [Ux_restr,Uy_restr]
Nodes[6].set_Restr([1,1])
for i in range(len(Elements)):
    N1=Elements[i][0]
    N2=Elements[i][1]
    Elem[i+1]=Bar_element(Nodes[N1],Nodes[N2],E=2e8,A=0.005)


#Step 4: Node numbering with structural dofs, forming association matrices
Count_P=0
DOF_counter=0
for i in Nodes.keys():
    #Take node
    node :Node=Nodes[i]
    restr=node.get_restrains()
    #Check its restraints , restr[j]==0 implies unrestrained, restr[j]==1 implies restrained
    #restr=node.get_Restr()
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


#Step 5: ASSEMBLING
    #Global stiffness maTriX [K_total] (2n x 2n), [Force Vector](nx1) and [Displacement Vector](nx1)
N=2*len(Nodes)
P=Count_P            #No. of Unrestrained joints for KPP :P x P
K_Ts=np.zeros((N,N)) #Total stiffness matrix
for i in Elem.keys():
    member :Bar_element=Elem[i]
    Association=member.get_Association()
    k_g=member.get_global_k()
    K_Ts[np.ix_(Association,Association)]=np.add(K_Ts[np.ix_(Association,Association)],k_g)

#print(K_Ts)

F=np.zeros((N,1))    #Nodal Force vector
U=np.zeros((N,1))    #Nodal Displacement Vector
for i in Nodes.keys():
    node :Node =Nodes[i]
    Association=node.get_Association()
    node_disp=node.get_Disp().reshape(-1,1)
    node_force=node.get_Load().reshape(-1,1)
    F[Association]=node_force
    U[Association]=node_disp

#Step 6: Partition [K_total] into Kpp,Kpx,Kxp,Kxx and [U] into Up,Ux, and [F] into Fp,Fx

Kpp=K_Ts[:P,:P]
Kpx=K_Ts[:P,P:]
Kxp=K_Ts[P:,:P]
Kxx=K_Ts[P:,P:]

Fp=F[:P]
Fx=F[P:]
Up=U[:P]
Ux=U[P:]

#Step 7: Solve the F-D relation equations
"""Solving force-displacement equations
    [Kpp]{Up}+[Kpx]{Ux}={Fp}
    [Kxp]{Up}+[Kxx]{Ux}={Fx}

    {Up}=[Kpp]^(-1){[Fp]-[Kpx][Ux]}     #Unknown forces 
    {Fx}=[Kxp]{Up}+[Kxx]{Ux}            #Unknown Reactions
"""
Up=np.dot(np.linalg.inv(Kpp),np.add(Fp,-1*np.dot(Kpx,Ux)))
Fx=np.add(np.dot(Kxp,Up),np.dot(Kxx,Ux))
U=np.concatenate((Up,Ux))
F=np.concatenate((Fp,Fx))

#Step 8: UPDATING NODAL VALUES with outputs
for i in Nodes.keys():
    node :Node =Nodes[i]
    Association=node.get_Association()
    node.set_Disp(U[Association])
    node.set_Load(F[Association])
    

#Step 9 :Print Results
print("Results")
print("\nNode\tUx(mm)\tUy(mm)\tFx(kN)\tFy(kN)")
for i in Nodes.keys():
    node :Node =Nodes[i]
    Disp_node=node.get_Disp()
    Disp_force=node.get_Load()
    print(i,"\t%8.7f\t %8.7f\t %8.7f\t %8.7f"%(Disp_node[0][0]*1000,Disp_node[1][0]*1000,Disp_force[0][0],Disp_force[1][0]))

print("\nElement\tStress(kN/m^2)")
for i in Elem.keys():
    member :Bar_element =Elem[i]
    stress=member.get_stress()
    Nodal_force=member.get_Nodal_Force()
    print(i,'\t',"%5.4f"%stress[0])


# %% [markdown]
# <H1>BONUS QUESTION</H1>
# Plot the deformed shape of truss obtained in Q10

# %%
#write solution here
Nod_cord=[(0,0),(5,7),(5,0),(10,7),(10,0),(15,0)]
New_cord=[]
Scale=2000
for i in Nodes.keys():
    coord=list(Nodes[i].get_new_pos(Scale).flat)
    New_cord+=[(coord[0],coord[1])]
Elements=[(1,2),(1,3),(2,3),(2,4),(2,5),(3,5),(4,5),(4,6),(5,6)]

fig,ax=plt.subplots(1,1)
ax.scatter([i[0] for i in New_cord],[i[1] for i in New_cord],s=100,color='k')

for i in Elements:
    ax.plot([New_cord[i[0]-1][0],New_cord[i[1]-1][0]],[New_cord[i[0]-1][1],New_cord[i[1]-1][1]],color='k',linestyle='--')
    ax.plot([Nod_cord[i[0]-1][0],Nod_cord[i[1]-1][0]],[Nod_cord[i[0]-1][1],Nod_cord[i[1]-1][1]],color='b')

ax.legend(["Deflected","Original"])
fig.savefig("2018CE10169_Q11.png")
plt.show()


# %%



