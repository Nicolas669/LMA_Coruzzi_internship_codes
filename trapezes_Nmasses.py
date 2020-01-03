_author_ = "Nicolas Coruzzi"
_filename_ = "trapezes_Nmasses"
_creationdate_ = "13/12/19"

import numpy as np

def trapezes_Nmasses(M,Utz,Vtz,t_tot,nb,N):

    Id=np.identity(2*N)
    tps=[]
    delta_t=t_tot/nb
    A1=Id-((delta_t*M)/2)
    A2=Id+((delta_t*M)/2)
    A1inverse=np.linalg.inv(A1)

    solus = np.zeros((2 * N, nb + 1), float)
    solus[(N - 1)][0] = Utz
    solus[(2 * N - 1)][0] = Vtz

    tn=0
    tps+=[tn]
    for n in range(nb):
        tn = tn + delta_t
        tps+=[tn]
        solus[:, n+1] = np.dot(np.dot(A1inverse,A2),solus[:, n])
    return (solus,tps)