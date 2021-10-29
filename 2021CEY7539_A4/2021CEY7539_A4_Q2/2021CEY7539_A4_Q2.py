#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display,Math


# # Renumber nodes 1,2,3 as 2,1,3

# In[2]:


def FrameElemStiffness(E,A,H,I):
    return np.array([[E*A/H,0,0,-E*A/H,0,0],[0,12*E*I/(H**3),6*E*I/(H**2),0,-12*E*I/(H**3),6*E*I/(H**2)],[0,6*E*I/(H**2),4*E*I/H,0,-6*E*I/(H**2),2*E*I/(H)],[-E*A/H,0,0,E*A/H,0,0],[0,-12*E*I/(H**3),-6*E*I/(H**2),0,12*E*I/(H**3),-6*E*I/(H**2)],[0,6*E*I/(H**2),2*E*I/(H),0,-6*E*I/(H**2),4*E*I/(H)]])


# In[3]:


display(Math('All\space distances\space are\space in\space mm\space and\space loads\space in\space kN.'))
display(Math('Areas\space are\space in\space mm^{2}\space and\space so\space on.'))
L1,L2,L3,w,E,A,I=input("Enter the values of L1,L2,L3,w,E (in GPa),A and I  in the same order separated by commas\n").split(',')
L1,L2,L3,w,E,A,I=float(L1),float(L2),float(L3),float(w),float(E),float(A),float(I)


# # Element stiffness matrices

# In[4]:


#Element 1                        (connecting nodes 2 and 1)
L_3=np.sqrt(L1**2+L3**2)
c1,s1=-L1/L_3,L3/L_3
k1_l=FrameElemStiffness(E,A,L_3,I)
T1=np.array([[c1,s1,0,0,0,0],[-s1,c1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,c1,s1,0],[0,0,0,-s1,c1,0],[0,0,0,0,0,1]])

k1=np.linalg.multi_dot([np.matrix.transpose(T1),k1_l,T1])

As_1=[1,2,3,4,5,6]

k2=FrameElemStiffness(E,A,L2,I)

As_2=[1,2,3,7,8,9]


# # Total structure stiffness matrix

# In[6]:


As_com=[As_1,As_2]
k_com=[k1,k2]

Kts=np.zeros([9,9])
for i in range(2):       # for 2 elements
    for j in range(6):   # both 6s for member dofs          
        for k in range(6):
            Kts[As_com[i][j]-1][As_com[i][k]-1]+=k_com[i][j][k]

#Correcting for renumbering
r=[4,5,6,1,2,3,7,8,9]

Kts_actual=np.zeros([9,9])

for i in range(9):
    for j in range(9):
        Kts_actual[r[i]-1][r[j]-1]+=Kts[i][j]
display(Math('K_{TS}\space(in\space kN/mm)=%s'%(sym.latex(sym.sympify(Kts_actual)))))


# # Kpp, Kxp and F

# In[11]:


Kpp=Kts[0:3,0:3]
Kxp=Kts[3:9,0:3]

F=np.transpose(np.array([[0,-w*L2/2,-(w*L2**2)/12]]))


# # Displacements

# In[12]:


Ikpp=np.linalg.inv(Kpp)                     

Dp=np.dot(Ikpp,F)
display(Math('Displacements\space at\space node\space 2 (in\space mm\space and\space rad)=%s'%sym.latex(sym.sympify(Dp))))


# # Reactions

# In[13]:


X_n=np.dot(Kxp,Dp)                    # This time, X includes forces acting along reactions
X=X_n+np.transpose(np.array([[0,0,0,0,w*L2/2,-(w*L2**2)/12]]))

display(Math('Reactions\space (in\space kN\space and\space kNmm)=%s'%sym.latex(sym.sympify(X))))


# # Element forces

# In[18]:


D1=np.zeros([6,1])
for i in range(6):
    if As_1[i]<=3:
        D1[i]+=Dp[As_1[i]-1]

D1_l=np.matmul(T1,D1)                                   # l for local

#local force vector
f_1=np.matmul(k1_l,D1_l)

D2=np.zeros([6,1])
for i in range(6):
    if As_1[i]<=3:
        D2[i]+=Dp[As_2[i]-1]

# Considering contribution from fixed-end forces

f_2=np.matmul(k2,D2)+np.transpose(np.array([[0,w*L2/2,(w*L2**2)/12,0,w*L2/2,-(w*L2**2)/12]]))
display(Math('Axial\space force(kN),\space shear\space force(kN)\space and\space bending\space moment(kNmm)\space for\space element\space 1 =%s'%sym.latex(sym.sympify(f_1))))
display(Math('Axial\space force(kN),\space shear\space force(kN)\space and\space bending\space moment(kNmm)\space for\space element\space 2 =%s'%sym.latex(sym.sympify(f_2))))


# # AFD (considering compression as negative)

# In[20]:


plt.subplot(1,2,1)
plt.plot([0,L_3],[-f_1[0],-f_1[0]])
plt.title('AFD for element 1')

plt.subplot(1,2,2)
plt.plot([0,L2],[-f_2[0],-f_2[0]])
plt.title('AFD for element 2')


plt.tight_layout()
plt.show()


# # SFD

# In[22]:


plt.subplot(1,2,1)
plt.plot([0,L_3],[f_1[1],f_1[1]])
plt.title('SFD for element 1')

plt.subplot(1,2,2)
plt.plot([0,L2],[f_2[1],f_2[1]])
plt.title('SFD for element 2')


plt.tight_layout()
plt.show()


# # BMD

# In[23]:


x1_plot=np.linspace(0,L_3,20)
x2_plot=np.linspace(0,L2,20)


M1_plot=-f_1[2]+f_1[1]*x1_plot
M2_plot=-f_2[2]+f_2[1]*x2_plot-w*x2_plot*x2_plot/2


plt.subplot(1,3,1)
plt.plot(x1_plot,M1_plot,'r-')
plt.title('BMD for element 1')

plt.subplot(1,3,2)
plt.plot(x2_plot,M2_plot,'k-')
plt.title('BMD for element 2')



plt.tight_layout()
plt.show()

