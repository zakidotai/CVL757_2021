#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import sympy as sym
from IPython.display import Math, display
import matplotlib.pyplot as plt


# In[3]:


display(Math('All\space distances\space are\space in\space mm\space and\space loads\space in\space kN.'))
display(Math('Areas\space are\space in\space mm^{2}\space and\space so\space on.'))
L1,L2,L3,w,E,I=input("Enter the values of L1, L2,L3, w (in kN/mm), E (in GPa) and I  in the same order separated by commas\n").split(',')
L1,L2,L3,w,E,I=float(L1),float(L2),float(L3),float(w),float(E),float(I)


# # Member stiffness matrices

# In[4]:


# Member 1 (connecting nodes 1 and 2)            

k1=np.array([[12*E*I/(L1**3),6*E*I/(L1**2),-12*E*I/(L1**3),6*E*I/(L1**2)],[6*E*I/(L1**2),4*E*I/L1,-6*E*I/(L1**2),2*E*I/(L1)],[-12*E*I/(L1**3),-6*E*I/(L1**2),12*E*I/(L1**3),-6*E*I/(L1**2)],[6*E*I/(L1**2),2*E*I/(L1),-6*E*I/(L1**2),4*E*I/(L1)]])

        
#Association matrix

As_1=[4,1,5,2]               #3 comes after 4 to account for the numbering of rotation first at the right end

# Member 2 (containing nodes 2 and 3)

k2=np.array([[12*E*I/(L2**3),6*E*I/(L2**2),-12*E*I/(L2**3),6*E*I/(L2**2)],[6*E*I/(L2**2),4*E*I/L2,-6*E*I/(L2**2),2*E*I/(L2)],[-12*E*I/(L2**3),-6*E*I/(L2**2),12*E*I/(L2**3),-6*E*I/(L2**2)],[6*E*I/(L2**2),2*E*I/(L2),-6*E*I/(L2**2),4*E*I/(L2)]])

As_2=[5,2,6,3]

#Member 3

k3=np.array([[12*E*I/(L3**3),6*E*I/(L3**2),-12*E*I/(L3**3),6*E*I/(L3**2)],[6*E*I/(L3**2),4*E*I/L3,-6*E*I/(L3**2),2*E*I/(L3)],[-12*E*I/(L3**3),-6*E*I/(L3**2),12*E*I/(L3**3),-6*E*I/(L3**2)],[6*E*I/(L3**2),2*E*I/(L3),-6*E*I/(L3**2),4*E*I/(L3)]])

As_3=[6,3,7,8]


# # Total structure stiffness matrix

# In[5]:


As_com=[As_1,As_2,As_3]
k_com=[k1,k2,k3]

Kts=np.zeros([8,8])
for i in range(3):       # for 3 elements
    for j in range(4):   # both 4s for member dofs          
        for k in range(4):
            Kts[As_com[i][j]-1][As_com[i][k]-1]+=k_com[i][j][k] 

display(Math('K_{TS}\space(in\space kN/mm)=%s'%(sym.latex(sym.sympify(Kts))))) 


# # Defining F, Kpp and Kxp

# In[6]:


Kpp=Kts[0:3,0:3]
Kxp=Kts[3:8,0:3]

# convert the UDL into equivalent joint loads
F=np.transpose(np.array([[0,(-w*L2**2)/12,(w*L2**2)/12]]))


# # Finding rotations

# In[7]:


Ikpp=np.linalg.inv(Kpp)                     #Inverse of Kpp as F=Kpp*Dp
Dp=np.dot(Ikpp,F)                           #Displacements
display(Math('Rotations\space at\space joints\space 1,2,\space and\space 3\space (in\space rad)=%s'%sym.latex(sym.sympify(Dp))))
print(f'Hence, vertical displacement and rotation at node 2 are {Dp[0]}mm and {Dp[1]}rad')
print(f'Rotation at node 3= {Dp[2]}rad')


# # Finding reactions

# In[8]:


X_n=np.dot(Kxp,Dp)                    # This time, X includes forces acting along reactions
X=X_n+np.transpose(np.array([[0,w*L2/2,w*L2/2,0,0]]))
display(Math('Reactions\space (in\space kN\space and\space kNmm)=%s'%sym.latex(sym.sympify(X))))


# # Shear forces in elements

# In[9]:


#Element 1
V1=X[0]
plt.subplot(1,3,1)
plt.plot([0,L1],[V1,V1])
plt.title('SFD for element 1')
plt.xlabel('$x(mm)$',fontsize=15)
plt.ylabel('$V(kN)$',fontsize=15)

#Element 2
x2_plot=np.linspace(0,L2,20)
V2=V1-w*x2_plot+X[1]
plt.subplot(1,3,2)
plt.plot(x2_plot,V2)
plt.title('SFD for element 2')
plt.xlabel('$x(mm)$',fontsize=15)
plt.ylabel('$V(kN)$',fontsize=15)

#Element 3
V3=X[3]
plt.subplot(1,3,3)
plt.plot([0,L3],[V3,V3])
plt.title('SFD for element 3')
plt.xlabel('$x(mm)$',fontsize=15)
plt.ylabel('$V(kN)$',fontsize=15)

plt.tight_layout()
plt.show()


# # Bending moments in elements

# In[10]:


x1_plot=np.linspace(0,L1,20)
x3_plot=np.linspace(0,L3,20)          

#Element 1
M1_plot=V1*x1_plot

#Element 2
M2_plot=V1*(L1+x2_plot)+(X[1]*x2_plot)-(w*x2_plot*x2_plot/2)

#Element 3
M3_plot=X[3]*(L3-x3_plot)+X[4]

plt.subplot(1,3,1)
plt.plot(x1_plot,M1_plot,'k-')
plt.xlabel('$x(mm)$', fontsize=15)
plt.ylabel('$M(kNmm)$',fontsize=15)
plt.title('BMD for element 1')
plt.subplot(1,3,2)
plt.plot(x2_plot,M2_plot,'b-')
plt.xlabel('$x(mm)$', fontsize=15)
plt.ylabel('$M(kNmm)$',fontsize=15)
plt.title('BMD for element 2')
plt.subplot(1,3,3)
plt.plot(x3_plot,M3_plot,'g-')
plt.xlabel('$x(mm)$', fontsize=15)
plt.ylabel('$M(kNmm)$',fontsize=15)
plt.title('BMD for element 3')
plt.tight_layout()
plt.show()

