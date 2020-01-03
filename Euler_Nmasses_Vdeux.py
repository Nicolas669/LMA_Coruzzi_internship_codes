_author_ = "Nicolas Coruzzi"
_filename_ = "Euler_Nmasses_Vdeux"
_creationdate_ = "13/12/19"

import numpy as np

def euler_exp_Nmasses(M,Utz,Vtz,t_tot,nb,N):

    tps=[]
    delta_t=t_tot/nb

    solus = np.zeros((2 * N, nb + 1), float)
    solus[(N - 1)][0] = Utz
    solus[(2 * N - 1)][0] = Vtz
    tn = 0
    tps += [tn]
    for n in range(nb):
        tn = tn + delta_t
        tps+=[tn]
        solus[:, n+1] = solus[:, n] + delta_t*(np.dot(M,solus[:,n]))
    return (solus,tps)

def euler_imp_Nmasses(M,Utz,Vtz,t_tot,nb,N):

    tps = []
    delta_t = t_tot / nb

    Id=np.identity(2*N)
    A=Id-delta_t*M
    Ainverse=np.linalg.inv(A)

    solus = np.zeros((2 * N, nb + 1), float)
    solus[(N - 1)][0] = Utz
    solus[(2 * N - 1)][0] = Vtz
    tn=0
    tps+=[tn]
    for n in range(nb):
        tn = tn + delta_t
        tps+=[tn]
        solus[:, n+1] = np.dot(Ainverse,solus[:, n])
    return (solus,tps)