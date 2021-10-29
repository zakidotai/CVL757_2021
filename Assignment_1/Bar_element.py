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
    