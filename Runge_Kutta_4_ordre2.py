_author_ = "Nicolas Coruzzi"
_filename_ = "Runge_Kutta_4_ordre2"
_creationdate_ = "02/12/19"

import numpy as np

def runge_kutta_ordre2_homogene(M,Xz,vz,t_tot,nb):

    Id=np.identity(2,float)
    tps=[]
    delta_t=t_tot/nb
    solus = np.zeros((2, nb+1),float)
    solus[:, 0] = (Xz,vz)
    tn=0
    tps+=[tn]

    M2=np.dot(M,M)
    M3=np.dot(M2,M)

    mat2=Id+ (delta_t*M)/2
    mat3=Id+(delta_t*M)/2+ ((delta_t/2)**2)*M2
    mat4=Id+ delta_t*M+ ((delta_t**2)/2)*M2+((delta_t**3)/4)*M3

    for n in range(nb):
        tn = tn + delta_t
        tps+=[tn]
        solus[:, n+1] = solus[:, n] + (delta_t/6)*(np.dot(M,solus[:,n])+2*np.dot(np.dot(M,mat2),solus[:,n])+2*np.dot(np.dot(M,mat3),solus[:,n])+np.dot(np.dot(M,mat4),solus[:,n]))
    return (solus,tps)
