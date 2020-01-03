_author_ = "Nicolas Coruzzi"
_filename_ = "trapezes_ordre2_homogene"
_creationdate_ = "01/12/19"

from fonctions_matrice import*

def trapezes_ordre2_homogene(M,Xz,vz,t_tot,nb):

    Id=np.identity(2)
    tps=[]
    delta_t=t_tot/nb
    A1=Id-((delta_t*M)/2)
    A2=Id+((delta_t*M)/2)
    A1inverse=inverse(A1)
    solus = np.zeros((2, nb+1))
    solus[:, 0] = (Xz,vz)
    tn=0
    tps+=[tn]
    for n in range(nb):
        tn = tn + delta_t
        tps+=[tn]
        solus[:, n+1] = np.dot(np.dot(A1inverse,A2),solus[:, n])
    return (solus,tps)