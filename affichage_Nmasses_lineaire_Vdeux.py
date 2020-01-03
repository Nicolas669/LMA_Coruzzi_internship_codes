_author_ = "Nicolas Coruzzi"
_filename_ = "affichage_Nmasses_lineaire_Vdeux"
_creationdate_ = "13/12/19"

#N=n+1 masses m
#ressorts de raideur k

#chaine de masses-ressorts (en lineaire)

from Euler_Nmasses_Vdeux import*
from trapezes_Nmasses import*
from Runge_Kutta_4_Nmasses import*
from Newmark_Nmasses import*
import matplotlib.pyplot as plt
import numpy as np


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
#mil[0]=1.

fin=np.array([-1]*(N-1),float)
#fin[-1]=0

nu = tridiag(np.array([-1]*(N-1),float), mil, fin)

periode=(2*np.pi)/((raid / m) ** (1 / 2))

nb_periode_t_tot=6
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

#############Proprietes de la fenetre
#if N==1:
#    plt.figure(figsize=(20,12), dpi=80)
    #plt.suptitle('Comparaison des methodes analytique et numerique sur une ED homogene du second ordre',fontsize=30)

    #plt.plot(temps_exp,sol_exp[:N][0],".",fillstyle='none',label='euler_exp',
    #         color="green")
#    plt.plot(temps_imp,sol_imp[:N][0],"+",fillstyle='none',label='euler_imp',
#             color="blue")
#    plt.plot(temps_trap,sol_trap[:N][0],".",fillstyle='none',label='trapezes',
#             color="black")
#    plt.plot(temps_rk,sol_rk[:N][0],"x",fillstyle='none',label='runge_kutta',
#             color="deeppink")
#    plt.plot(temps_nmark,sol_nmark[0],"v", fillstyle='none',label='newmark',
#             color="goldenrod")

#    plt.legend(bbox_to_anchor=(0.54, 0.7), loc='lower left',
#               fontsize =14,borderaxespad=0.1)
#    plt.ylabel('position de la masse',fontsize=16)
#    plt.xlabel('temps (secondes)',fontsize=16)
    #plt.xscale('log')
#    plt.savefig("test_1masse")
#    plt.show()
#if N==2:

#    plt.figure(figsize=(12, 8), dpi=80)
#    plt.plot(temps_imp, sol_imp[:N][:N][0], ".", fillstyle='none', label='euler_imp_masse1',
#             color="blue")
#    plt.plot(temps_imp, sol_imp[:N][:N][1], "x", fillstyle='none', label='euler_imp_masse2',
#             color="blue")
#    plt.legend(bbox_to_anchor=(0.54, 0.7), loc='lower left',
#               fontsize=14, borderaxespad=0.1)
#    plt.ylabel('position des masses', fontsize=16)
#    plt.xlabel('temps (secondes)', fontsize=16)
#    plt.savefig("euler_imp_2masses")

#    plt.figure(figsize=(12, 8), dpi=80)
#    plt.plot(temps_trap, sol_trap[:N][:N][0], ".", fillstyle='none', label='trapezes_masse1',
#             color="black")
#    plt.plot(temps_trap, sol_trap[:N][:N][1], "x", fillstyle='none', label='trapezes_masse2',
#             color="black")
#    plt.legend(bbox_to_anchor=(0.54, 0.7), loc='lower left',
#               fontsize=14, borderaxespad=0.1)
#    plt.ylabel('position des masses', fontsize=16)
#    plt.xlabel('temps (secondes)', fontsize=16)
#    plt.savefig("trap_2masses")

#    plt.figure(figsize=(12, 8), dpi=80)
#    plt.plot(temps_rk, sol_rk[:N][:N][0], ".", fillstyle='none', label='runge_kutta_masse1',
#             color="deeppink")
#    plt.plot(temps_rk, sol_rk[:N][:N][1], "x", fillstyle='none', label='runge_kutta_masse2',
#             color="deeppink")
#    plt.legend(bbox_to_anchor=(0.54, 0.7), loc='lower left',
#               fontsize=14, borderaxespad=0.1)
#    plt.ylabel('position des masses', fontsize=16)
#    plt.xlabel('temps (secondes)', fontsize=16)
#    plt.savefig("rk_2masses")

#    plt.figure(figsize=(12, 8), dpi=80)
#    plt.plot(temps_nmark, sol_nmark[0], ".", fillstyle='none', label='newmark_masse1',
#             color="goldenrod")
#    plt.plot(temps_nmark, sol_nmark[1], "x", fillstyle='none', label='newmark_masse2',
#             color="goldenrod")
#    plt.legend(bbox_to_anchor=(0.54, 0.7), loc='lower left',
#               fontsize=14, borderaxespad=0.1)
#    plt.ylabel('position des masses', fontsize=16)
#    plt.xlabel('temps (secondes)', fontsize=16)
#    plt.savefig("nmark_2masses")

#    plt.show()


plt.figure(figsize=(12, 8), dpi=80)
for p in range(N-1):
    plt.plot(temps_rk, sol_rk[:N][:N][p]+p*Utz+1, "x", fillstyle='none',
             color="deeppink")
plt.plot(temps_rk, sol_rk[:N][:N][N-1] + (N-1) * Utz+1, "x", fillstyle='none',label='Runge Kutta',
         color="deeppink")
plt.legend(bbox_to_anchor=(0.45, 0.92), loc='lower left',
               fontsize=14, borderaxespad=0.1)
plt.yticks(np.arange(1, N+1, step=Utz))
plt.ylabel('mass number', fontsize=16)
plt.xlabel('time (seconds)', fontsize=16)
plt.savefig("rk_masses"+str(N))

plt.figure(figsize=(12, 8), dpi=80)
for p in range(N-1):
    plt.plot(temps_trap, sol_trap[:N][:N][p]+p*Utz+1, ".", fillstyle='none',
             color="black")
plt.plot(temps_trap, sol_trap[:N][:N][N-1] + (N-1) * Utz+1, ".", fillstyle='none',label='Trapezoidal',
         color="black")
plt.legend(bbox_to_anchor=(0.45, 0.92), loc='lower left',
               fontsize=14, borderaxespad=0.1)
plt.yticks(np.arange(1, N+1, step=Utz))
plt.ylabel('mass number', fontsize=16)
plt.xlabel('time (seconds)', fontsize=16)
plt.savefig("trap_masses"+str(N))

plt.figure(figsize=(12, 8), dpi=80)
for p in range(N-1):
    plt.plot(temps_imp, sol_imp[:N][:N][p]+p*Utz+1, "+", fillstyle='none',
             color="blue")
plt.plot(temps_imp, sol_imp[:N][:N][N-1] + (N-1) * Utz+1, "+", fillstyle='none',label='Euler implicite',
         color="blue")
plt.legend(bbox_to_anchor=(0.45, 0.92), loc='lower left',
               fontsize=14, borderaxespad=0.1)
plt.yticks(np.arange(1, N+1, step=Utz))
plt.ylabel('mass number', fontsize=16)
plt.xlabel('time (seconds)', fontsize=16)
plt.savefig("imp_masses"+str(N))

plt.figure(figsize=(12, 8), dpi=80)
for p in range(N-1):
    plt.plot(temps_nmark, sol_nmark[p]+p*Utz+1, "v", fillstyle='none',
             color="goldenrod")
plt.plot(temps_nmark, sol_nmark[N-1] + (N-1) * Utz+1, "v", fillstyle='none',label='Newmark',
         color="goldenrod")
plt.legend(bbox_to_anchor=(0.45, 0.92), loc='lower left',
               fontsize=14, borderaxespad=0.1)
plt.yticks(np.arange(1, N+1, step=Utz))
plt.ylabel('mass number', fontsize=16)
plt.xlabel('time (seconds)', fontsize=16)
plt.savefig("nmark_masses"+str(N))

plt.show()
