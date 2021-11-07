# Q1 : Reversing List
l = [1, 2, 3, 4, 5, 6, 7, 8]
Input  = [12, 35, 9, 56, 24]
def reverse_list(l):
    return l[::-1]
print(f"Original List is : {l}")
print(f"Reverse List is : {reverse_list(l)}")
print(f"Original List is : {Input}")
print(f"Reverse List is : {reverse_list(Input)}")

# Q2 : Swapping Elements in Lists
l =[1, 2, 3, 4, 5, 6, 7, 8]
p1 = 3
p2 = 5
def swapFunc(l,p1,p2):
    l[p1], l[p2] = l[p2], l[p1]
    return l
print(f"Original List is : {l}")
print(f"Swapped List is : {swapFunc(l,(p1-1),(p2-1))}") #**SInce in indexing of list starts from 0 we make p1-1 for correction**

 
# Q3: Reversing a string
s = "CVL756"
def reverse_string(s):
    return s[::-1]
print(f"Original String is : {s}")
print(f"Reverse String is : {reverse_string(s)}")


# Q4: Reversing the order of words in a string

s = "This is CVL756"
name= 'Pratik Kamble'
def reverse_sentence(s):
    rev = (' '.join(reversed(s.split(' ')))) 
    # in above sentence we first split sentence in words and then reverse and join with <space>
    return rev
print(f"Original Sentence is : {s}")
print(f"Reverse sentence is : {reverse_sentence(s)}")
print(f"Original Sentence is : {name}")
print(f"Reverse sentence is : {reverse_sentence(name)}")

# Q5: Initialise a 3x3 matrix with random values, use numpy
import numpy as np
np.random.seed(20)
random_matrix = np.random.rand(3,3)
print(random_matrix)

# Q6: Swap the columns of random_matrix
random_matrix_col = np.array([[0.5881308,  0.89771373, 0.89153073],
                              [0.81583748, 0.03588959, 0.69175758],
                              [0.37868094, 0.51851095, 0.65795147]])
c1 = 1
c2 = 2
def swap_column(l,c1,c2):
    l[:, [c1, c2]] = l[:, [c2, c1]]
    return l
print(f"Matrix with swapped column :\n {swap_column(random_matrix_col,c1-1,c2-1)}")

# Q7: Swap the rows of random_matrix
random_matrix_row = np.array([[0.5881308,  0.89771373, 0.89153073],
                              [0.81583748, 0.03588959, 0.69175758],
                              [0.37868094, 0.51851095, 0.65795147]])
r1 = 1
r2 = 2
def swap_rows(l,r1,r2):
    l[[r1]], l[[r2]] = l[[r2]], l[[r1]]
    return l
print(f"Matrix with swapped row :\n {swap_rows(random_matrix_row,r1-1,r2-1)}")


# Q8: Plot the sine and cosine curves as shown below
import matplotlib.pyplot as plt
import numpy as np
x = np.arange(0,2*np.pi,0.001)
y1 = np.sin(x)
y2 = np.cos(x)
plt.plot(x,y1,'b', label='sine(\u03B8)')
plt.plot(x,y2,'r',label='cosine(\u03B8)')
plt.xlabel('\u03B8 values from 0 to 2\u03C0')
plt.ylabel('sin(\u03B8)\n cos(\u03B8)')
plt.grid()
plt.legend(loc= "best")
plt.show()

# Q9: Plot the truss
import numpy as np
import matplotlib.pyplot as plt
ndcon= np.array([[1,2],[1,3],[3, 2],[2, 4],[2, 5],[3, 5],[5, 4],[4, 6],[5,6]])
nelem = len(ndcon)
x=[0,5,5,10,10,15]
y=[0,7,0,7,0,0]
ielem = 0
while ielem < 9: 
    x1 =x[ndcon[ielem][0]-1]
    x2 =x[ndcon[ielem][1]-1]
    y1 =y[ndcon[ielem][0]-1]
    y2 =y[ndcon[ielem][1]-1]
    plt.plot([x1,x2],[y1,y2],'Black',marker='o',markerfacecolor='r',markersize='12')
    ielem+=1
plt.show()  

# Q.10 Solution
import numpy as np
A= 0.005  # meter^2
E = 2e8   # KN/m^2
L1=5      # meter
L2=7      # meter
L3=8.6023 # meter
ndcon= np.array([[1,2],[1,3],[3, 2],[2, 4],[2, 5],[3, 5],[5, 4],[4, 6],[5,6]])
nelem = len(ndcon)
x=[0,5,5,10,10,15]
y=[0,7,0,7,0,0]
nodes = len(x)
ndofn=2
nnode = 2
tdofs = nodes*ndofn
gstiff = np.zeros([tdofs,tdofs])
import math as mt
ielem=0
while ielem < 9:   
    x1 =x[ndcon[ielem][0]-1]
    x2 =x[ndcon[ielem][1]-1]
    y1 =y[ndcon[ielem][0]-1]
    y2 =y[ndcon[ielem][1]-1]
    gbdof=[]
    L1 = mt.sqrt((x2-x1)**2+(y2-y1)**2)
    C =(x2-x1)/L1
    S = (y2-y1)/L1
    A1 = A*E/L1
    K1 =A1* np.matrix('%s, %s , %s, %s; %s, %s, %s, %s; %s, %s, %s, %s; %s, %s, %s, %s' % (C**2, C*S,-C**2, -C*S, C*S, S**2, -C*S, -S**2,-C**2, -C*S, C**2, C*S, -C*S, -S**2, C*S, S**2))
    
    inode =0
    while inode < nnode:
        idofn=1
        while idofn <= ndofn:
            gbdof.append((ndcon[ielem,inode]-1)*ndofn+idofn)
            idofn += 1
        inode += 1
    i = 0
    if i == 0:
        gstiff[gbdof[i]-1,gbdof[i]-1]=gstiff[gbdof[i]-1,gbdof[i]-1] + K1[0,0]
        gstiff[gbdof[i]-1,gbdof[i+1]-1]=gstiff[gbdof[i]-1,gbdof[i+1]-1] + K1[0,1]
        gstiff[gbdof[i]-1,gbdof[i+2]-1]=gstiff[gbdof[i]-1,gbdof[i+2]-1] + K1[0,2]
        gstiff[gbdof[i]-1,gbdof[i+3]-1]=gstiff[gbdof[i]-1,gbdof[i+3]-1] + K1[0,3]
    i = 1
    if i == 1:
        gstiff[gbdof[i]-1,gbdof[i-1]-1]=gstiff[gbdof[i]-1,gbdof[i-1]-1] + K1[1,0]
        gstiff[gbdof[i]-1,gbdof[i]-1]=gstiff[gbdof[i]-1,gbdof[i]-1]     + K1[1,1]
        gstiff[gbdof[i]-1,gbdof[i+1]-1]=gstiff[gbdof[i]-1,gbdof[i+1]-1] + K1[1,2]
        gstiff[gbdof[i]-1,gbdof[i+2]-1]=gstiff[gbdof[i]-1,gbdof[i+2]-1] + K1[1,3]
    i = 2
    if i == 2:
        gstiff[gbdof[i]-1,gbdof[i-2]-1]=gstiff[gbdof[i]-1,gbdof[i-2]-1] + K1[2,0]
        gstiff[gbdof[i]-1,gbdof[i-1]-1]=gstiff[gbdof[i]-1,gbdof[i-1]-1] + K1[2,1]
        gstiff[gbdof[i]-1,gbdof[i]-1]=gstiff[gbdof[i]-1,gbdof[i]-1]     + K1[2,2]
        gstiff[gbdof[i]-1,gbdof[i+1]-1]=gstiff[gbdof[i]-1,gbdof[i+1]-1] + K1[2,3]
    i = 3 
    if i == 3:
        gstiff[gbdof[i]-1,gbdof[i-3]-1]=gstiff[gbdof[i]-1,gbdof[i-3]-1] + K1[3,0]
        gstiff[gbdof[i]-1,gbdof[i-2]-1]=gstiff[gbdof[i]-1,gbdof[i-2]-1] + K1[3,1]
        gstiff[gbdof[i]-1,gbdof[i-1]-1]=gstiff[gbdof[i]-1,gbdof[i-1]-1] + K1[3,2]
        gstiff[gbdof[i]-1,gbdof[i]-1]=gstiff[gbdof[i]-1,gbdof[i]-1]     + K1[3,3]
    ielem += 1
    
print("Global Stiffness Matrix = ")
print(gstiff)

b = np.delete(gstiff,0, 0)
b = np.delete(b,0,1)
b = np.delete(b,0, 0)
b = np.delete(b,0,1)
b = np.delete(b,8, 0)
b = np.delete(b,8,1)
b = np.delete(b,8, 0)
b = np.delete(b,8,1)
force=np.matrix('0;0;20;0;0;0;0;0;0;0;0;0')
f1 = np.delete(force,0, 0)
f1 = np.delete(f1,0,0)
f1= np.delete(f1,8, 0)
f1= np.delete(f1,8,0)
u = np.linalg.inv(b).dot(f1)
U=np.zeros([tdofs,1])

i=0
while i<8:
    U[i+2][0]=u[i][0]
    i +=1
print("Displacement at Joints")    
print(np.matrix([U[0],U[1],U[10],U[11]]))
print("other displacement")
print(u)
ielem=0
print("Positive sign = Tensile Force \nNegative Sign = Compression")
while ielem < 9: 
    x1 =x[ndcon[ielem][0]-1]
    x2 =x[ndcon[ielem][1]-1]
    y1 =y[ndcon[ielem][0]-1]
    y2 =y[ndcon[ielem][1]-1] 
    L1 = mt.sqrt((x2-x1)**2+(y2-y1)**2)
    C =(x2-x1)/L1
    S = (y2-y1)/L1
    stress=(E/L1)*np.matrix([-C,-S,C,S]).dot(np.matrix([U[ndcon[ielem][0]*2-2],U[ndcon[ielem][0]*2-1],U[ndcon[ielem][1]*2-2],U[ndcon[ielem][1]*2-1]]))
    print(f"stress in element {ielem+1} is : {stress} KN/m^2")
    ielem +=1
print("\n Calculation of reactions at node 1 and 6")
# Calculation of the reactions    
F = np.matrix(gstiff).dot(np.matrix(U))
print(f"Horizontal Reaction at node 1 = {F[0]}\nVertical Reactiona at node 1 = {F[1]}\nHorizontal Reaction at node 6 = {F[10]}\nVertical Reactiona at node 6 = {F[11]}")


# Plot the defomred shape of truss obtained in Q10
#**Dispacement values are too small to visible in the plot so that to plot the deformed shape the displacement are magnified**
import numpy as np
import matplotlib.pyplot as plt
ndcon= np.array([[1,2],[1,3],[3, 2],[2, 4],[2, 5],[3, 5],[5, 4],[4, 6],[5,6]])
nelem = len(ndcon)
x=[0,5+2.18759959e-0,5+ 1.11111111e-01,10+ 1.85426626e-0,10+ 2.22222222e-01,15]
y=[0,7+-3.50052910e-01,0-3.50052910e-01,7+ 1.11957672e-01,0+-5.41375661e-01,0]
x1=[0,5,5,10,10,15]
y1=[0,7,0,7,0,0]
def TrussPlot(x,y,line_style):
    ielem = 0
    while ielem < 9: 
        x1 =x[ndcon[ielem][0]-1]
        x2 =x[ndcon[ielem][1]-1]
        y1 =y[ndcon[ielem][0]-1]
        y2 =y[ndcon[ielem][1]-1]
        line_style=str(line_style)
        plt.plot([x1,x2],[y1,y2],line_style)
        ielem+=1
TrussPlot(x,y,'r')
TrussPlot(x1,y1,'b')
plt.show()