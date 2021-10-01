#!/usr/bin/env python
# coding: utf-8

# <h1>Assignment 1

# <h3>Q1: Reversing a list

# In[1]:


#following are sample input and corresponding outputs
Input  = [12, 35, 9, 56, 24]
Output = [24, 35, 9, 56, 12]

Input  = [1, 2, 3]
Output = [3, 2, 1]


# In[4]:


#write your solution here to reverse a list
I=[12,35,9,56,24]
I.reverse()
print(I)
I=[1,2,3]
I.reverse()
print(I)


# <h3>Q2: Swapping elements in lists

# In[7]:


Input  = [12, 35, 9, 56, 24]
Output = [12, 9, 35, 56, 24]


# In[8]:


#write your solution here
temp=Input[1]
Input[1]=Input[2]
Input[2]=temp
print(Input)


# <h3>Q3: Reversing a string

# In[ ]:


input_str  = 'This is CVL757'
output_str = '757LVC si sihT'


# In[12]:


#write your solution here
s='This is CVL757'
st=""
for i in s:
    st=i+st
print (st)


# <h3>Q4: Reversing the order of words in a string

# In[13]:


input_str  = 'This is CVL757'
output_str = 'CVL757 is This'


# In[19]:


#write your solution here
sentence='This is CVL757'
sentbreak=sentence.split(' ')
print(sentbreak)
revSent=""
for i in sentbreak:
    revSent= i+' '+revSent
print(revSent)


# <h3>Q5: Initialise a 3x3 matrix with random values, use numpy

# In[23]:


#random_matrix = #write solution here
import numpy as np
A=np.random.random((3,3))
print (A)


# In[ ]:





# <h3>Q6: Swap the columns of random_matrix

# In[35]:


inp = np.array([[0.13806091, 0.82474243, 0.15646752],
                [0.50029824, 0.45875794, 0.53401557],
                [0.95397773, 0.79407795, 0.07442586]])

out = np.array([[0.15646752, 0.82474243, 0.13806091],
                [0.53401557, 0.45875794, 0.50029824],
                [0.07442586, 0.79407795, 0.95397773]])


# In[36]:


#write your solution here to swap 1st and last column
temp = np.copy(inp[:,0])
inp[:,0]=inp[:,2]
inp[:,2]=temp
print(inp)


# <h3>Q7: Swap the rows of random_matrix

# In[37]:


inp = np.array([[0.13806091, 0.82474243, 0.15646752],
               [0.50029824, 0.45875794, 0.53401557],
               [0.95397773, 0.79407795, 0.07442586]])

out = np.array([[0.95397773, 0.79407795, 0.07442586],
               [0.50029824, 0.45875794, 0.53401557],
               [0.13806091, 0.82474243, 0.15646752]])


# In[38]:


#write your solution here to swap 1st and 3rd rows
temp=np.copy(inp[0,:])
inp[0,:]=inp[2,:]
inp[2,:]=temp
print (inp)


# <h3>Q8: Plot the sine and cosine curves as shown below

# In[40]:


import matplotlib.pyplot as plt


# In[41]:


x = np.arange(0,2*np.pi,0.01)
#write code here to plot and show the legend also
y=np.sin(x)
z=np.cos(x)
plt.plot(x,y,x,z)
plt.legend(['sin(x)','cos(x)'])
plt.show()


# <h3>Q9: Plot the following truss</h3>
# Length of element 2, 6 and 9 (between nodes 1 and 3, 3 and 5, and 5 and 9) is 5.<br>
# Length of element 3 and 7 is 7m.
#  

# <img src='q9.svg' alt="Truss Q9" width="500" height="600">

# In[54]:


#write solution here
points = np.array([[0,0],[5,7],[5,0], [10,7], [10,0], [15,0]])
lines =np.array([[points[0], points[1]],
                 [points[0], points[2]],
                 [points[1], points[2]],
                 [points[1], points[3]],
                 [points[1], points[4]],
                 [points[2], points[5]],
                 [points[3], points[4]],
                 [points[3], points[5]],
                 [points[4], points[5]]])
for line in lines:
    plt.plot(line[:,0], line[:,1])
plt.show()


# **Q10: Consider the plane truss shown above. Given E = 200GPa and A = 0.005m2, and horizontal load of 20kN at node 2. Both node 1 and node 6 have pin supports:**
# 
# 1. the global stiffness matrix for the structure.
# 2. the horizontal and vertical displacements at nodes 2, 3, 4, and 5.
# 3. the horizontal and vertical reactions at nodes 1 and 6.
# 4. the stress in each element.

# In[5]:


#write solution here
import numpy as np
nodes=[[0,0],
                 [5,7],
                 [5,0],
                 [10,7],
                 [10,0],
                 [15,0]]
members=             [[1,2],
                     [1,3],
                     [3,2],
                     [2,4],
                     [5,2],
                     [3,5],
                     [5,4],
                     [6,4],
                     [5,6]]
dof=2*len(nodes)
structure_stiffness=np.zeros((dof,dof))
for i in range(0,9):
    
    fnode=members[i][0]
    snode=members[i][1]
    fnodecordinate=members[fnode-1][:]
    snodecordinate=members[snode-1][:]
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


# <H1>BONUS QUESTION</H1>
# Plot the deformed shape of truss obtained in Q10

# In[7]:


#write solution here
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
    plt.plot(x,y,color='black',marker='o', markerfacecolor='black', markersize=1)
    
    x=[fnodecordinate[0],snodecordinate[0]]
    y=[fnodecordinate[1],snodecordinate[1]]
    plt.plot(x,y,color='red',linestyle='dashed',marker='o', markerfacecolor='red', markersize=1)


# In[ ]:





# In[ ]:





# In[ ]:




