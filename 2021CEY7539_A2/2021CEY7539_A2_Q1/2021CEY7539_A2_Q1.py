#!/usr/bin/env python
# coding: utf-8

# In[75]:


import numpy as np
from IPython.display import display, Math
import sympy as sym
import matplotlib.pyplot as plt


# # Renumbering the nodes
# ## Nodes 1,2,3 are renumbered as 3,1,2
# ### This allows us the partition the total structure stiffness matrix in a desired manner
# #### Degrees of freedom at each node are numbered in the following order: Vertical, anticlockwise rotation except at the hinge support where the rotation is numbered first

# In[27]:


display(Math('All\space distances\space are\space in\space mm\space and\space loads\space in\space kN.'))
display(Math('Areas\space are\space in\space mm^{2}\space and\space so\space on.'))
L1,L2,P,E,I=input("Enter the values of L1, L2, P, E (in GPa) and I  in the same order separated by commas\n").split(',')
L1,L2,P,E,I=float(L1),float(L2),float(P),float(E),float(I)


# # Member stiffness matrices

# In[37]:


# Member 1 (connecting nodes 1 and 2)            #renumbered nodes

k1=np.array([[12*E*I/(L2**3),6*E*I/(L2**2),-12*E*I/(L2**3),6*E*I/(L2**2)],[6*E*I/(L2**2),4*E*I/L2,-6*E*I/(L2**2),2*E*I/(L2)],[-12*E*I/(L2**3),-6*E*I/(L2**2),12*E*I/(L2**3),-6*E*I/(L2**2)],[6*E*I/(L2**2),2*E*I/(L2),-6*E*I/(L2**2),4*E*I/(L2)]])
        
#Association matrix

As_1=[1,2,4,3]               #3 comes after 4 to account for the numbering of rotation first at the right end

# Member 2 (containing nodes 3 and 1)

k2=np.array([[12*E*I/(L1**3),6*E*I/(L1**2),-12*E*I/(L1**3),6*E*I/(L1**2)],[6*E*I/(L1**2),4*E*I/L1,-6*E*I/(L1**2),2*E*I/(L1)],[-12*E*I/(L1**3),-6*E*I/(L1**2),12*E*I/(L1**3),-6*E*I/(L1**2)],[6*E*I/(L1**2),2*E*I/(L1),-6*E*I/(L1**2),4*E*I/(L1)]])

As_2=[5,6,1,2]
display(Math('k_{1}=k_{2}\space (in\space kN/mm)=%s  '%sym.latex(sym.sympify(k1))))


# # Total structure stiffness matrix (KTS)

# In[39]:


k_com=[k1,k2]

As_com=[As_1,As_2]

Kts=np.zeros([6,6])
for i in range(2):       # for 2 elements
    for j in range(4):   # both 4s for member dofs          
        for k in range(4):
            Kts[As_com[i][j]-1][As_com[i][k]-1]+=k_com[i][j][k] 

#Correcting for reordering   nodes 1,2,3= actualnodes 2,3,1

r=[3,4,5,6,1,2]

Kts_actual=np.zeros([6,6])

for i in range(6):
    for j in range(6):
        Kts_actual[r[i]-1][r[j]-1]+=Kts[i][j]

display(Math('K_{TS}\space(in\space kN/mm)=%s'%(sym.latex(sym.sympify(Kts_actual)))))         


# # Defining force vector F, Kpp and Kxp

# In[43]:


Kpp=Kts[0:3,0:3]
Kxp=Kts[3:6,0:3]
F=np.transpose(np.array([[-P,0,0]]))


# # Finding displacements and rotations

# In[73]:


Ikpp=np.linalg.inv(Kpp)                     #Inverse of Kpp as F=Kpp*Dp
Dp=np.dot(Ikpp,F)                           #Displacements
display(Math('Displacements\space (in\space mm\space and\space rad)=%s'%sym.latex(sym.sympify(Dp))))
print(f'Hence, vertical displacement and rotation at node 2 are {Dp[0]}mm and {Dp[1]}rad')
print(f'Rotation at node 3= {Dp[2]}rad')


# # Finding reactions

# In[72]:


X=np.dot(Kxp,Dp)                    # This time, X includes forces acting along reactions
display(Math('Reactions\space (in\space kN\space and\space kNmm)=%s'%sym.latex(sym.sympify(X))))
print(f'Reactions at node 1 are {X[1]}kN and{X[2]}kNmm')
print(f'Reaction at node 3 = {X[0]}kN')


# # Shear forces in elements

# In[147]:


# Element 1
V1=-X[0]
V2=X[1]
print(f'Shear force in right element={V1}kN\nShear force in left element={V2}kN')
plt.subplot(1,2,2)
plt.plot([0,1],[V1,V1])
plt.axis([0,1,2*V1,0])
plt.title('SFD for the right element')
plt.tick_params(bottom=False,labelbottom=False)
ax1=plt.gca()
ax1.spines['right'].set_color('none')              #removes the right spine(boundary line of the figure)
# ax1.spines['top'].set_color('none')
# ax1.spines['left'].set_color('none')              
ax1.spines['bottom'].set_color('none')
plt.subplot(1,2,1)
plt.plot([0,1],[V2,V2])
plt.axis([0,1,0,2*V2])
ax2=plt.gca()
plt.title('SFD for the left element')
plt.tick_params(bottom=False,labelbottom=False)
ax2.spines['right'].set_color('none')              #removes the right spine(boundary line of the figure)
ax2.spines['top'].set_color('none')
# ax2.spines['left'].set_color('none')              
# ax2.spines['bottom'].set_color('none')
plt.show()


# # Bending moments in elements

# In[146]:


x1_plot=np.linspace(0,L2,20)
x2_plot=np.linspace(0,L1,20)
#Element 1

M1_plot=V2*(L1+x1_plot)-X[2]-P*x1_plot

#Element 2
M2_plot=(V2*x2_plot)-X[2]

plt.subplot(1,2,2)
plt.plot(x1_plot,M1_plot,'k-')
plt.xlabel('$x(mm)$', fontsize=15)
plt.ylabel('$M_{right}(kNmm)$',fontsize=15)
plt.title('BMD for right element')
plt.subplot(1,2,1)
plt.plot(x2_plot,M2_plot,'b-')
plt.xlabel('$x(mm)$', fontsize=15)
plt.ylabel('$M_{left}(kNmm)$',fontsize=15)
plt.title('BMD for left element')
plt.tight_layout()
plt.show()


# In[144]:


[5]-np.array([5,6,8,-8])


# In[ ]:




