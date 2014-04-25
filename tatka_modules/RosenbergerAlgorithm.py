# -*- coding:utf-8-*
"""Algorithm S detection"""
import numpy as np

#Define I
I = np.identity(3)


#How to go from u(n) to u(n+1)
def update_(U, D, d, landa):
    m = U.T*d
    p = (I-U*U.T)*d
    p_norm = np.linalg.norm(p)
    
    U_left = np.hstack((U,1/p_norm*p))
    Q = np.hstack((landa*D,m))
    Q = np.vstack((Q,np.matrix([0,0,p_norm])))

    #SVD
    U_right, D_new, V_left = np.linalg.svd(Q)

    #Get rid of the smaller eigen value
    D_new = D_new[0:2]
    D_new = np.diagflat(D_new)

    U_right = U_right[:,0:2]

    return U_left*U_right, D_new


#Algorithm that return S and P wave separatly
def rosenberger(dataX, dataY, dataZ, landa):
    
    #Construct the matrix
    A = np.matrix(np.vstack((dataZ,dataX,dataY)))

    # SVD of the first 3 samples:
    U, D, V = np.linalg.svd(A[:,0:3])

    # Get rid of the smallest eigen value
    D = D[0:2]
    D = np.diagflat(D)
    U = U[:,0:2]

    save_U = np.zeros(len(dataX)-2)
    save_U[0] = abs(U[0,0])

    #Initialysing two matrix (what is the assumption here?)
    Dp = np.matrix(np.zeros((3, len(dataX)-2)))
    Ds = np.matrix(np.zeros((3, len(dataX)-2)))

    Dp[:,0] = abs(U[0,0])*A[:,2]
    Ds[:,0] = (1-abs(U[0,0]))*A[:,2]
 
    #loop over all the values
    for i in range(0,A.shape[1]-3):
        d = A[:,i+3]
        U, D = update_(U, D, d, landa)
  
        Dp[:,i+1] = abs(U[0,0])*d
        Ds[:,i+1] = (1-abs(U[0,0]))*d
    
        save_U[i+1] = abs(U[0,0])

    return np.array(Dp), np.array(Ds), save_U
