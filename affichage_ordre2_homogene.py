_author_ = "Nicolas Coruzzi"
_filename_ = "affichage_ordre2_homogene"
_creationdate_ = "01/12/19"

from analytique_ordre2_homogene import*
from Euler_ordre2_homogene import*
from trapezes_ordre2_homogene import*
from Runge_Kutta_4_ordre2 import*
from Newmark_ordre2_homogene import*
import matplotlib.pyplot as plt

#raid=float(input("Saisir la raideur du ressort (N/m) : "))
#long=float(input("Saisir la longueur a vide du ressort (m) : "))
#m=float(input("Saisir la masse (kg) : "))
#ca=float(input("Saisir la valeur du coefficient d'amortissement (kg/s) : "))
#Xz=float(input("Saisir la valeur de X au temps zero (m) : "))
#vz=float(input("Saisir la valeur de X prime au temps zero (m/s): "))
#t_tot=int(input("Saisir la duree (entier) (s) : "))
#nb=int(input("Saisir le nombre de points (entier positif): "))

raid=5
long=1
m=1
ca=0.3
Xz=1
vz=0

periode=(2*np.pi)/((raid / m) ** (1 / 2))

nb_periode_t_tot=5
t_tot=nb_periode_t_tot*periode
nb=50*nb_periode_t_tot

omega = (raid / m) ** (1 / 2)
beta = ca / (m * 2 * omega)

M = np.array([[0, 1], [-(omega ** 2), -2 * omega * beta]])

#t_tot doit être quelques fois la periode (avec periode= 2*np.pi/omega)
#nb de sorte à avoir environ 30 pts par periodes

(sol_ana,temps_ana)=analytique_ordre2_homogene_complete(omega,beta,Xz,vz,t_tot,nb)
(sol_exp,temps_exp)=euler_exp_ordre2_homogene(M,Xz,vz,t_tot,nb)
(sol_imp,temps_imp)=euler_imp_ordre2_homogene(M,Xz,vz,t_tot,nb)
(sol_trap,temps_trap)=trapezes_ordre2_homogene(M,Xz,vz,t_tot,nb)
(sol_rk,temps_rk)=runge_kutta_ordre2_homogene(M,Xz,vz,t_tot,nb)
(sol_nmark,temps_nmark,vitesse_nmark)=newmark_ordre2_homogene(raid,m,ca,Xz,vz,t_tot,nb)

#############Proprietes de la fenetre

plt.figure(figsize=(20,12), dpi=80)
#plt.suptitle('Comparaison des methodes analytique et numerique sur une ED homogene du second ordre',fontsize=30)

#plt.annotate(' coefficient amortissement : '+ str(ca)+
#             ' \n position initiale : '+ str(Xz)+
#             ' \n vitesse initiale : '+ str(vz),
#             fontsize=20, xy = (0.66,0.78),
#             xycoords='figure fraction', xytext = (0.66,0.78),
#             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
#                           'width': 15, 'headwidth': 30},
#             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
#                       edgecolor="forestgreen", lw=1,))


plt.plot(temps_ana,sol_ana,"-",label='analytic',
         color="firebrick")
#plt.plot(temps_exp,sol_exp[0],".",fillstyle='none',label='euler_exp',
#         color="green")
plt.plot(temps_imp,sol_imp[0],"+",fillstyle='none',label='euler implicite',
         color="blue")
plt.plot(temps_trap,sol_trap[0],".",fillstyle='none',label='trapezoidal',
         color="black")
plt.plot(temps_rk,sol_rk[0],"x",fillstyle='none',label='Runge Kutta',
         color="deeppink")
plt.plot(temps_nmark,sol_nmark,"v",fillstyle='none',label='Newmark',
         color="goldenrod")
plt.legend(bbox_to_anchor=(0.35, 0.8), loc='lower left',
           fontsize =14,borderaxespad=0.1)

plt.ylabel('mass position',fontsize=16)
plt.xlabel('time (seconds)',fontsize=16)
#plt.xscale('log')
plt.savefig("homogene_ordre2")
plt.show()