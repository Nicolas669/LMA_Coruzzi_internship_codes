_author_ = "Nicolas Coruzzi"
_filename_ = "affichage_1masse_nonlineaire"
_creationdate_ = "18/12/19"

import matplotlib.pyplot as plt
from Euler_1masse_nonlineaire import*
from trapezes_1masse_nonlineaire import*
from runge_kutta_1masse_nonlineaire import*


raid=5.
m=1.
ca=0.35
xz=1.
nb=500
t_tot=10
vz=0.


#pour le non lineaire:
#xz=5.
#nb=3000
#t_tot=2


#pour newton raphson
precision=10**(-10)
nb_iteration=10
####################


def fonction(x,k):
    """renvoie la fonction et sa dérivee"""
    fon=k*(x+x**3)
    der=(1+3*(x**2))*k
    return (fon,der)

def FONCTION(elt,alph,m,fonc,k):
    F1 = elt[1]
    F2 = ((-alph / m) * elt[1]) - ((fonc(elt[0],k)[0]) / m)
    F = np.array([F1, F2])
    JacF = np.array([[0., 1.], [-((fonc(elt[0],k)[1]) / m), (-alph / m)]])
    return (F, JacF)

#(sol_ana,temps_ana)=analytique_ordre2_homogene_complete(omega,beta,Xz,vz,t_tot,nb)
(sol_imp,temps_imp)=euler_imp_1masse_nonlineaire(xz,vz,t_tot,nb,FONCTION,ca,m,fonction,raid,precision,nb_iteration)
(sol_trap,temps_trap)=trapeze_1masse_nonlineaire(xz,vz,t_tot,nb,FONCTION,ca,m,fonction,raid,precision,nb_iteration)
(sol_rk,temps_rk)=runge_kutta_1masse_nonlineaire(xz,vz,t_tot,nb,FONCTION,ca,m,fonction,raid,precision,nb_iteration)

#############Proprietes de la fenetre

plt.figure(figsize=(20,12), dpi=80)
#plt.suptitle('Comparaison des methodes analytique et numerique sur une ED homogene du second ordre',fontsize=30)

plt.plot(temps_imp,sol_imp[0],"+",fillstyle='none',label='Euler implicite',
         color="blue")
plt.plot(temps_trap,sol_trap[0],".",fillstyle='none',label='trapezoidal',
         color="black")
plt.plot(temps_rk,sol_rk[0],"x",fillstyle='none',label='Runge Kutta',
         color="deeppink")

plt.legend(bbox_to_anchor=(0.46, 0.7), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('mass position',fontsize=16)
plt.xlabel('time (seconds)',fontsize=16)
#plt.xscale('log')
plt.savefig("1masse_nonlineaire")
plt.show()