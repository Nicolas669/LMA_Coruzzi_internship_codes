_author_ = "Nicolas Coruzzi"
_filename_ = "Euler_ordre2_homogene"
_creationdate_ = "01/12/19"

from fonctions_matrice import*

def euler_exp_ordre2_homogene(M,Xz,vz,t_tot,nb):
    tps=[]
    delta_t=t_tot/nb
    #if delta_t >= (2 /max_mod_val_prop(M)):
    #    print("probleme, voir euler_exp_ordre2_homogene")
    solus = np.zeros((2, nb+1),float)
    solus[:, 0] = (Xz,vz)
    tn=0
    tps+=[tn]
    for n in range(nb):
        tn = tn + delta_t
        tps+=[tn]
        solus[:, n+1] = solus[:, n] + delta_t*(np.dot(M,solus[:,n]))
    return (solus,tps)

def euler_imp_ordre2_homogene(M,Xz,vz,t_tot,nb):

    Id=np.identity(2)
    tps=[]
    delta_t=t_tot/nb
    A=Id-delta_t*M
    Ainverse=inverse(A)
    solus = np.zeros((2, nb+1))
    solus[:, 0] = (Xz,vz)
    tn=0
    tps+=[tn]
    for n in range(nb):
        tn = tn + delta_t
        tps+=[tn]
        solus[:, n+1] = np.dot(Ainverse,solus[:, n])
    return (solus,tps)

