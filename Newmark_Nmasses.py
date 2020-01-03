_author_ = "Nicolas Coruzzi"
_filename_ = "Newmark_Nmasses"
_creationdate_ = "13/12/19"

import numpy as np


def newmark_Nmasses(omega,beta,nu, Utz, Vtz, t_tot, nb, N):

    #dans ce programme, D est la position, V la vitesse et A l'acceleration

    tps = []
    delta_t = t_tot / nb

    M=np.identity(N,float)

    C=2*beta*omega*np.identity(N,float)

    K=(omega**2)*nu

    Gamma=1/2
    Beta=1/4

    solus = np.zeros((2 * N, nb + 1), float)
    solus[(N - 1)][0] = Utz
    solus[(2 * N - 1)][0] = Vtz
    D = solus[:N]
    V = solus[N:2*N]

    A=np.zeros((N,nb+1),float)
    A[:,0]=-np.dot(C,V[:,0])-np.dot(K,D[:,0])

    S = M + Gamma*delta_t*C + Beta*(delta_t**2)*K

    tn=0
    tps+=[tn]
    for n in range(1,nb+1):
        tn = tn + delta_t
        tps+=[tn]

        Dpr= D[:,(n-1)]+delta_t*V[:,(n-1)]+(delta_t**2)*((1/2)-Beta)*A[:,n-1]
        Vpr= V[:,(n-1)]+delta_t*(1-Gamma)*A[:,n-1]

        theta= -np.dot(C,Vpr) -np.dot(K,Dpr)

        A[:,n]=np.dot(theta,np.linalg.inv(S))

        D[:,n]=Dpr+Beta*(delta_t**2)*A[:,n]
        V[:,n]=Vpr+Gamma*delta_t*A[:,n]
    return (D,tps,V)