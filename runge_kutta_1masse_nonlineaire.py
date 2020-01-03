_author_ = "Nicolas Coruzzi"
_filename_ = "runge_kutta_1masse_nonlineaire"
_creationdate_ = "19/12/19"

from Newton_Raphson import*

def runge_kutta_1masse_nonlineaire(xz,vz,t_tot,nb,FONC,alph,m,fonc,k,precision,nb_iteration):
    """alpha et le coef d'amortissement.
        fonc est notre fonction non lineaire, qui renvoie la valeur de
        notre fonction [0] et de sa derivee [1] en un point"""

    tps=[]
    delta_t=t_tot/nb
    Y = np.zeros((2, nb+1))
    Y[:, 0] = (xz,vz)

    def G_rk(U,U_n,FONC):
        U_n1 = U_n
        U_n2 = U_n1 + (delta_t /2)*(FONC(U_n1,alph,m,fonc,k))[0]
        U_n3 = U_n1 + (delta_t /2)*(FONC(U_n2,alph,m,fonc,k))[0]
        U_n4 = U_n1 + delta_t * (FONC(U_n3, alph, m, fonc, k))[0]

        G=U-U_n-(delta_t/6)*(FONC(U_n1,alph,m,fonc,k)[0]+2*FONC(U_n2,alph,m,fonc,k)[0]+2*FONC(U_n3,alph,m,fonc,k)[0]+FONC(U_n4,alph,m,fonc,k)[0])
        dG=np.identity(2,float)-(delta_t/6)*(FONC(U_n1,alph,m,fonc,k)[1]+2*FONC(U_n2,alph,m,fonc,k)[1]+2*FONC(U_n3,alph,m,fonc,k)[1]+FONC(U_n4,alph,m,fonc,k)[1])
        return (G,dG)

    tn=0
    tps+=[tn]
    for n in range(nb):
        tn = tn + delta_t
        tps+=[tn]
        Y[:, n+1] = newton_raphson_vectoriel(precision,nb_iteration,Y[:,0],G_rk,Y[:,n],FONC)
    return (Y,tps)