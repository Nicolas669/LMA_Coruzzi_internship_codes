_author_ = "Nicolas Coruzzi"
_filename_ = "Euler_1masse_nonlineaire"
_creationdate_ = "18/12/19"

from Newton_Raphson import*

def euler_imp_1masse_nonlineaire(xz,vz,t_tot,nb,FONC,alph,m,fonc,k,precision,nb_iteration):
    """alpha et le coef d'amortissement.
        fonc est notre fonction non lineaire, qui renvoie la valeur de
        notre fonction [0] et de sa derivee [1] en un point"""

    tps=[]
    delta_t=t_tot/nb
    Y = np.zeros((2, nb+1))
    Y[:, 0] = (xz,vz)

    def G_euler_imp(U,U_n,FONC):
        G=U-delta_t*(FONC(U,alph,m,fonc,k))[0]-U_n
        dG=np.identity(2,float)-delta_t*(FONC(U,alph,m,fonc,k))[1]
        return (G,dG)

    tn=0
    tps+=[tn]
    for n in range(nb):
        tn = tn + delta_t
        tps+=[tn]
        Y[:, n+1] = newton_raphson_vectoriel(precision,nb_iteration,Y[:,0],G_euler_imp,Y[:,n],FONC)
    return (Y,tps)
