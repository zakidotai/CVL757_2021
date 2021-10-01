#!/usr/bin/env python
# coding: utf-8

# In[5]:


def Reverse(lst):
   lst.reverse()
   return lst
    
lst = [12, 35, 9, 56, 24]
lst1=[1,2,3]
print(Reverse(lst))
print(Reverse(lst1))

input=[12,35,9,56,24]
b=input[1]
input[1]=input[2]
input[2]=b
print(input)

input= 'This is CVL757'
n=len(input)
#print(n)
for i in range(n):
 print(input[n-1-i], end =" ")

input = 'This is CVL757'
word=input.split()
word.reverse()
print(' '.join(word))

# np.random.seed(100)
# random_matrix=np.random.rand(3,3)
# print (random_matrix)

import numpy as np
np.random.seed(100)
random_matrix=np.random.rand(3,3)
print (random_matrix)
random_matrix[:, [2, 0]] = random_matrix[:, [0, 2]]
print (random_matrix)

import numpy as np
np.random.seed(100)
random_matrix=np.random.rand(3,3)
print (random_matrix)
random_matrix[[2, 0],: ] = random_matrix[[0, 2], :]
print (random_matrix)

import matplotlib.pyplot as plt
import numpy as np
x= np.arange(0,2*np.pi,0.01)
plt.plot(x,np.sin(x),label="sin(x)")
plt.plot(x,np.cos(x),label="cos(x)")
plt.title("sin(x) and cos(x)")
plt.legend()
plt.xlabel("angle (x) in radians")
plt.ylabel("sin(x) and cos(x)")                         
plt.show()

import matplotlib.pyplot as plt
import numpy as np
cordinate_matrix=[[0,0],
                 [5,7],
                 [5,0],
                 [10,7],
                 [10,0],
                 [15,0]]
#print(cordinate_matrix)
element_connectivity=[[1,2],
                     [1,3],
                     [3,2],
                     [2,4],
                     [5,2],
                     [3,5],
                     [5,4],
                     [6,4],
                     [5,6]]
#print(element_connectivity)

for i in range(0,len(element_connectivity)):
    fnode=element_connectivity[i][0]
    snode=element_connectivity[i][1]
    fnodecordinate=cordinate_matrix[fnode-1][:]
    snodecordinate=cordinate_matrix[snode-1][:]
    x=[fnodecordinate[0],snodecordinate[0]]
    y=[fnodecordinate[1],snodecordinate[1]]
    plt.plot(x,y,color='black',marker='o', markerfacecolor='black', markersize=5)
    plt.annotate('%s' % i,x,y)
plt.show()
#Q10   
import numpy as np

cordinate_matrix=[[0,0],
                 [5,7],
                 [5,0],
                 [10,7],
                 [10,0],
                 [15,0]]

element_connectivity=[[1,2],
                     [1,3],
                     [3,2],
                     [2,4],
                     [5,2],
                     [3,5],
                     [5,4],
                     [6,4],
                     [5,6]]

dof=2*len(cordinate_matrix)
structure_stiffness=np.zeros((dof,dof))

for i in range(0,len(element_connectivity)):
    
    fnode=element_connectivity[i][0]
    snode=element_connectivity[i][1]
    fnodecordinate=cordinate_matrix[fnode-1][:]
    snodecordinate=cordinate_matrix[snode-1][:]
    x1=fnodecordinate[0]
    y1=fnodecordinate[1]
    x2=snodecordinate[0]
    y2=snodecordinate[1]    
    length= (((x2-x1)**2)+((y2-y1)**2))**0.5
    #print('/n',length)
    c=(x2-x1)/length
    s=(y2-y1)/length
    if x1==x2:
      tan=0
    else:
      tan=s/c
      
    Transformation_matrix=[[c, s, 0,0],[0,0,c,s]]
    Transformation_matrix = np.array(Transformation_matrix)
    transpose_Transformation_matrix=np.transpose(Transformation_matrix)
    k=0.005*2*(10**11)/length
    local_element_stiffness=[[k,-k],[-k,k]]
    local_element_stiffness = np.array(local_element_stiffness)
    global_element_stiffness=transpose_Transformation_matrix@local_element_stiffness@Transformation_matrix
    #print('local_element_stiffness',local_element_stiffness)
    #print('global_element_stiffness',global_element_stiffness)
   
    expanded_element_stiffness=np.zeros((dof,dof))
                
    expanded_element_stiffness[2*fnode-2,2*fnode-2]=global_element_stiffness[0,0]
    expanded_element_stiffness[2*fnode-2,2*fnode-1]=global_element_stiffness[0,1]
    expanded_element_stiffness[2*fnode-2,2*snode-2]=global_element_stiffness[0,2]
    expanded_element_stiffness[2*fnode-2,2*snode-1]=global_element_stiffness[0,3]
            
    expanded_element_stiffness[2*fnode-1,2*fnode-2]=global_element_stiffness[1,0]
    expanded_element_stiffness[2*fnode-1,2*fnode-1]=global_element_stiffness[1,1]
    expanded_element_stiffness[2*fnode-1,2*snode-2]=global_element_stiffness[1,2]
    expanded_element_stiffness[2*fnode-1,2*snode-1]=global_element_stiffness[1,3]
            
    expanded_element_stiffness[2*snode-2,2*fnode-2]=global_element_stiffness[2,0]
    expanded_element_stiffness[2*snode-2,2*fnode-1]=global_element_stiffness[2,1]
    expanded_element_stiffness[2*snode-2,2*snode-2]=global_element_stiffness[2,2]
    expanded_element_stiffness[2*snode-2,2*snode-1]=global_element_stiffness[2,3]
            
    expanded_element_stiffness[2*snode-1,2*fnode-2]=global_element_stiffness[3,0]
    expanded_element_stiffness[2*snode-1,2*fnode-1]=global_element_stiffness[3,1]
    expanded_element_stiffness[2*snode-1,2*snode-2]=global_element_stiffness[3,2]
    expanded_element_stiffness[2*snode-1,2*snode-1]=global_element_stiffness[3,3]
            
    structure_stiffness=structure_stiffness+expanded_element_stiffness

print('structure stiffness matrix: \n', structure_stiffness)
print('\n')

force=np.zeros(dof)
displacement=np.zeros(dof)
force[2]=20000
force1=force[2:10]
#print(force,force1)
stiffness=structure_stiffness[2:10,2:10]
#print("stiffness",stiffness)
inv_stiffness=np.linalg.inv(stiffness)
displacement1=inv_stiffness@force1
displacement[2:10]=displacement1                            
force=structure_stiffness@displacement

print('Reaction in Newton and Displacement in metre:\n')                  
print('rx1=',force[0])  
print('ry1=',force[1])  
print('rx6=',force[10])  
print('ry6=',force[11]) 
print('\n')
print('displacement of 2nd node in (x,y)=(',displacement[2],',',displacement[3],')')                   
print('displacement of 3nd node in (x,y)=(',displacement[4],',',displacement[5],')')   
print('displacement of 4th node in (x,y)=(',displacement[6],',',displacement[7],')')   
print('displacement of 5th node in (x,y)=(',displacement[8],',',displacement[9],')')  
print('\n')                            
strain=np.zeros(len(element_connectivity)) 
stress=np.zeros(len(element_connectivity)) 
global_element_displacement=np.zeros(4)
for i in range(0,9):
    fnode=element_connectivity[i][0]
    snode=element_connectivity[i][1]
    fnodecordinate=cordinate_matrix[fnode-1][:]                        
    snodecordinate=cordinate_matrix[snode-1][:]
    x1=fnodecordinate[0]
    y1=fnodecordinate[1]
    x2=snodecordinate[0]
    y2=snodecordinate[1]    
    length= (((x2-x1)**2)+((y2-y1)**2))**0.5
    c=(x2-x1)/length
    s=(y2-y1)/length
    Transformation_matrix=[[c, s, 0,0],
                           [0,0,c,s]]
    global_element_displacement[0]=displacement[2*fnode-2]
    global_element_displacement[1]=displacement[2*fnode-1]
    global_element_displacement[2]=displacement[2*snode-2]
    global_element_displacement[3]=displacement[2*snode-1]
    local_element_displacement=Transformation_matrix@global_element_displacement    
    strain[i]= (local_element_displacement[1]- local_element_displacement[0])*(1/length)
    stress[i]=2*(10**11)*strain[i]
    print('Stress(N/m^2) in members between node',fnode, 'and', snode, 'is',stress[i])  
#BONUS
import matplotlib.pyplot as plt
import numpy as np
cordinate_matrix=[[0,0],
                 [5,7],
                 [5,0],
                 [10,7],
                 [10,0],
                 [15,0]]

displacementxy=[[0,0],
                [displacement[2],displacement[3]],
                [displacement[4],displacement[5]],
                [displacement[6],displacement[7]],
                [displacement[8],displacement[9]],
                [0,0]]
displacementxy=np.multiply(displacementxy,1000)
print(displacementxy)
new_cordinate_matrix=np.add(cordinate_matrix,displacementxy)
print(new_cordinate_matrix)
print(cordinate_matrix)
element_connectivity=[[1,2],
                     [1,3],
                     [3,2],
                     [2,4],
                     [5,2],
                     [3,5],
                     [5,4],
                     [6,4],
                     [5,6]]

for i in range(0,len(element_connectivity)):
    fnode=element_connectivity[i][0]
    snode=element_connectivity[i][1]
    fnodecordinate=cordinate_matrix[fnode-1][:]
    snodecordinate=cordinate_matrix[snode-1][:]
    fnodecordinatenew=new_cordinate_matrix[fnode-1][:]
    snodecordinatenew=new_cordinate_matrix[snode-1][:]
    
    x=[fnodecordinatenew[0],snodecordinatenew[0]]
    y=[fnodecordinatenew[1],snodecordinatenew[1]]
    plt.plot(x,y,color='blue',marker='o', markerfacecolor='blue', markersize=1)
    
    x=[fnodecordinate[0],snodecordinate[0]]
    y=[fnodecordinate[1],snodecordinate[1]]
    plt.plot(x,y,color='black',linestyle='dashed',marker='o', markerfacecolor='black', markersize=1)
plt.show()


# In[ ]:




