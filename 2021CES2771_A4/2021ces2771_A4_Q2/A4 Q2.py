
import numpy as np
from scipy.linalg import block_diag
import matplotlib.pyplot as plt
E=200e9
I=1e-6
A=0.04
L1=2
L2=4
L3=3
w=8000
M=w*L2*L2/12
P=w*L2/2
KTS=np.zeros((9,9))
association=np.array([(1,2,3,4,5,6),(7,8,9,1,2,3)])
ff=np.array([(0,P,M,0,P,-M),(0,0,0,0,0,0)])
ff=np.transpose(ff)
P_mat=np.zeros(3)
P_mat[0]=0
P_mat[1]=-P
P_mat[2]=-M
P_mat=np.transpose(P_mat)
Fdash=np.array([0,-P,M,0,0,0])
Fdash=np.transpose(Fdash)
Ux=np.zeros(6)
Ux=np.transpose(Ux)
u=np.zeros(6)
Ke=np.zeros(((6,6,2)))
for i in range(2):
    if i==0:
        l=L2
        T=np.identity(6)
        Tt=np.transpose(T)
    else:
        l=(L1**2+L3**2)**0.5
        c=L1/l
        s=-L3/l
        t=s/c
        R=np.array([(c,s,0),(-s,c,0),(0,0,1)])
        T=block_diag(R,R)
        Tt=np.transpose(T)
    a=12*E*I/(l**3)
    b=6*E*I/(l**2)
    c=4*E*I/l
    d=2*E*I/l
    e=A*E/l
    Ke[:,:,i]=np.array([(e,0,0,-e,0,0),(0,a,b,0,-a,b),(0,b , c ,0, -b ,d),(-e,0,0,e,0,0),(0,-a, -b,0, a, -b),(0,b, d,0, -b, c)])
    Kg=Tt@Ke[:,:,i]@T
    for j in range(6):
        for k in range(6):
            KTS[association[i,j]-1,association[i,k]-1]= KTS[association[i,j]-1,association[i,k]-1]+Kg[j,k]
Kpp=KTS[0:3,0:3]
Kpx=KTS[0:3,3:9]
Kxp=KTS[3:9,0:3]
Kxx=KTS[3:9,3:9]
print('Global Stiffnes Matrix \n',KTS)
Kpp_inverse=np.linalg.inv(Kpp)
print(Kpp_inverse.shape)
print(Kpp_inverse.shape)
Up=Kpp_inverse@P_mat
X=Kxp@Up-Fdash
print('Fdash',Fdash.shape)
Displacement=np.concatenate((Up, Ux), axis=0)
Force=np.concatenate((P_mat, X), axis=0)

print('Force Matrix \n',Force)
print('Displacement Matrix \n',Displacement)
for i in range (2):
    for j in range(6):
        u[j]=Displacement[association[i,j]-1]
    if i==0:
        l=L2
        f=Ke[:,:,i]@u+ff[:,i]
        print('\n1st Member Force\n',f)
        print('\n1st Member Displacement\n',u)

        x1 = np.linspace(0, l,100)
        y1 = f[2]-f[1]*x1 +w* x1 * x1/2

        x2 = np.linspace(0, l,100)
        y2 = -f[1]+w*x2

        x3 = [0,l]
        y3 =  [-f[0],f[3]]

        fig, ((ax1, ax2, ax3)) = plt.subplots(nrows=1, ncols=3, sharex=True)
        fig.set_size_inches(14, 5)

        ax1.plot(x1, 0*x1, color='black', linewidth=2)
        ax1.plot(x1, y1, label='BMD', color="orangered", linewidth=2)  # Plot some data on the axes.
        ax1.set_xlabel('Length')  # Add an x-label to the axes.
        ax1.set_ylabel('BM in Nm')  # Add a y-label to the axes.
        ax1.set_title("BMD for ELEMENT 1")  # Add a title to the axes.
        ax1.legend(loc='upper left')  # Add a legend.
        ax1.grid(True)

        ax2.plot(x2, 0*x2, color='black', linewidth=2)
        ax2.plot(x2, y2, label='SFD', color="midnightblue", linewidth=2)  # Plot some data on the axes.
        ax2.set_xlabel('Length')  # Add an x-label to the axes.
        ax2.set_ylabel('shear force in N')  # Add a y-label to the axes.
        ax2.set_title("SFD for ELEMENT 1")  # Add a title to the axes.
        ax2.legend(loc='upper left')  # Add a legend.
        ax2.grid(True)

        ax3.plot(x3, [0,0], color='black', linewidth=2)
        ax3.plot(x3, y3, label='AFD', color="green", linewidth=2)  # Plot some data on the axes.
        ax3.set_xlabel('Length')  # Add an x-label to the axes.
        ax3.set_ylabel('Axial force in N')  # Add an x-label to the axes.
        ax3.set_title("AFD for ELEMENT 1")  # Add a title to the axes.
        ax3.legend(loc='upper right')  # Add a legend.
        ax3.grid(True)
        plt.show()

    else:
        l=(L1**2+L3**2)**0.5
        c=L1/l
        s=-L3/l
        t=s/c
        R=np.array([(c,s,0),(-s,c,0),(0,0,1)])
        T=block_diag(R,R)
        u=T@u
        f=Ke[:,:,i]@u+ff[:,i]
        print('\n2nd Member Force\n',f)
        print('\n2nd Member Displacement\n',u)
        x1= [0,l]
        y1= [f[2], f[5]]

        x2 = [0,l]
        y2 =  [f[1],f[4]]

        x3 = [0,l]
        y3 =  [-f[0],f[3]]

        fig, ((ax1, ax2, ax3)) = plt.subplots(nrows=1, ncols=3, sharex=True)
        fig.set_size_inches(14, 5)

        ax1.plot(x1, [0, 0], color='black', linewidth=2)
        ax1.plot(x1, y1, label='BMD', color="orangered", linewidth=2)  # Plot some data on the axes.
        ax1.set_xlabel('Length')  # Add an x-label to the axes.
        ax1.set_ylabel('BM in Nm')  # Add a y-label to the axes.
        ax1.set_title("BMD for ELEMENT 2")  # Add a title to the axes.
        ax1.legend(loc='upper left')  # Add a legend.
        ax1.grid(True)

        ax2.plot(x2, [0, 0], color='black', linewidth=2)
        ax2.plot(x2, y2, label='SFD', color="midnightblue", linewidth=2)  # Plot some data on the axes.
        ax2.set_xlabel('Length')  # Add an x-label to the axes.
        ax2.set_ylabel('shear force in N')  # Add a y-label to the axes.
        ax2.set_title("SFD for ELEMENT 2")  # Add a title to the axes.
        ax2.legend(loc='upper left')  # Add a legend.
        ax2.grid(True)

        ax3.plot(x3, [0, 0], color='black', linewidth=2)
        ax3.plot(x3, y3, label='AFD', color="green", linewidth=2)  # Plot some data on the axes.
        ax3.set_xlabel('Length')  # Add an x-label to the axes.
        ax3.set_ylabel('Axial force in N')  # Add an x-label to the axes.
        ax3.set_title("AFD for ELEMENT 2")  # Add a title to the axes.
        ax3.legend(loc='upper right')  # Add a legend.
        ax3.grid(True)
        plt.show()




