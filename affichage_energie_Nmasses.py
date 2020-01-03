_author_ = "Nicolas Coruzzi"
_filename_ = "affichage_energie_Nmasses"
_creationdate_ = "17/12/19"

from Euler_Nmasses_Vdeux import*
from trapezes_Nmasses import*
from Runge_Kutta_4_Nmasses import*
from Newmark_Nmasses import*
import matplotlib.pyplot as plt
import numpy as np

#Attention il faut ca=0 pour la conservation d'energie!
raid=5
long=1
m=1
ca=0.
Utz=1.
Vtz=0.
#Utz et Vtz sont ici la position et vitesse de la masse n+1 a l'instant initial (on suppose les autres au repos)

N=16
#N est le nombre de masse

#####################################################################
def tridiag(a, b, c):
    return np.diag(a, -1) + np.diag(b, 0) + np.diag(c, 1)

mil=np.array([2]*N,float)
mil[-1]=1.

fin=np.array([-1]*(N-1),float)

nu = tridiag(np.array([-1]*(N-1),float), mil, fin)

periode=(2*np.pi)/((raid / m) ** (1 / 2))

nb_periode_t_tot=5
t_tot=nb_periode_t_tot*periode
nb=50*nb_periode_t_tot

omega = (raid / m) ** (1 / 2)
beta = ca / (m * 2 * omega)

M1=np.concatenate((np.zeros([N,N],float), np.identity(N,float)), axis=1)
M2=np.concatenate((-(omega ** 2)*nu, (-2 * omega * beta)*np.identity(N,float)), axis=1)
M=np.concatenate((M1, M2), axis=0)

#####################################################################

(sol_imp,temps_imp)=euler_imp_Nmasses(M,Utz,Vtz,t_tot,nb,N)
(sol_trap,temps_trap)=trapezes_Nmasses(M,Utz,Vtz,t_tot,nb,N)
(sol_rk,temps_rk)=runge_kutta_Nmasses(M,Utz,Vtz,t_tot,nb,N)
(sol_nmark,temps_nmark,vitesse_nmark)=newmark_Nmasses(omega,beta,nu, Utz, Vtz, t_tot, nb, N)

#U1= sol_imp[:N][0]
#U1 point=solus2[N:2*N][0]
#U2= sol_imp[:N][1]
#U2 point=solus2[N:2*N][1]

energie_ana=[(1/2)*Vtz**2+(omega**2 / 2)*Utz**2]*len(temps_imp)

energie_imp=[]
energie_trap=[]
energie_rk=[]
energie_nmark=[]

for j in range(len(temps_imp)):
    r_imp = (omega ** 2 / 2) * sol_imp[:N][0][j] ** 2
    r_trap = (omega ** 2 / 2) * sol_trap[:N][0][j] ** 2
    r_rk = (omega ** 2 / 2) * sol_rk[:N][0][j] ** 2
    r_nmark=(omega ** 2 / 2) * sol_nmark[0][j] ** 2
    for k in range(N):
        r_imp += (1 / 2) * sol_imp[N:2 * N][k][j] ** 2
        r_trap += (1 / 2) * sol_trap[N:2 * N][k][j] ** 2
        r_rk += (1 / 2) * sol_rk[N:2 * N][k][j] ** 2
        r_nmark += (1 / 2) * vitesse_nmark[k][j] ** 2
    for l in range(N - 1):
        r_imp += (omega ** 2 / 2) * ((sol_imp[:N][l][j] - sol_imp[:N][l + 1][j]) ** 2)
        r_trap += (omega ** 2 / 2) * ((sol_trap[:N][l][j] - sol_trap[:N][l + 1][j]) ** 2)
        r_rk += (omega ** 2 / 2) * ((sol_rk[:N][l][j] - sol_rk[:N][l + 1][j])** 2)
        r_nmark+=(omega ** 2 / 2) * ((sol_nmark[l][j] - sol_nmark[l + 1][j])** 2)

    energie_imp+=[r_imp]
    energie_trap+=[r_trap]
    energie_rk+=[r_rk]
    energie_nmark += [r_nmark]

#############Proprietes de la fenetre

plt.figure(figsize=(12,8), dpi=80)
plt.plot(temps_imp,energie_ana,"-",label='analytic',
         color="firebrick")
plt.plot(temps_imp,energie_imp,"+",fillstyle='none',label='Euler implicite',
         color="blue")
plt.plot(temps_trap,energie_trap,".",fillstyle='none',label='trapezoidal',
         color="black")
plt.plot(temps_rk,energie_rk,"x",fillstyle='none',label='Runge Kutta',
         color="deeppink")
plt.plot(temps_nmark,energie_nmark,"v",fillstyle='none',label='Newmark',
         color="goldenrod")


plt.legend(bbox_to_anchor=(0.7, 0.5), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('System energy',fontsize=16)
plt.xlabel('time (seconds)',fontsize=16)
#plt.yscale('log')
plt.savefig("conservation_energie_"+str(N)+"masses")
plt.show()