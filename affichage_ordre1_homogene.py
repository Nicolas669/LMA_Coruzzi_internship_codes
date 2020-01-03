_author_ = "Nicolas Coruzzi"
_filename_ = "affichage_ordre1_homogene"
_creationdate_ = "01/12/19"

from analytique_ordre1_homogene import*
from Euler_ordre1_homogene import*
from trapezes_ordre1_homogene import*
from Runge_Kutta_4_ordre1 import*
import matplotlib.pyplot as plt
import numpy as np

c=2
xz=10

t_tot=5
nb=100

(sol_ana,temps_ana)=analytique_ordre1_homogene_complete(c,xz,t_tot,nb)
(sol_exp,temps_exp)=euler_exp_ordre1_homogene(c,xz,t_tot,nb)
(sol_imp,temps_imp)=euler_imp_ordre1_homogene(c,xz,t_tot,nb)
(sol_trap,temps_trap)=trapezes_ordre1_homogene(c,xz,t_tot,nb)
(sol_rk,temps_rk)=runge_kutta_ordre1_homogene(c,xz,t_tot,nb)

#############Proprietes de la fenetre

plt.figure(figsize=(20,12), dpi=80)
#plt.suptitle('Comparaison des methodes analytique et numeriques sur une ED homogene du premier ordre',fontsize=30)

plt.annotate(' coefficient c : '+ str(c)+
             ' \n position initiale : '+ str(xz)+
             ' \n nombre de points : '+ str(nb),
             fontsize=20, xy = (0.7,0.78),
             xycoords='figure fraction', xytext = (0.7,0.78),
             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
                           'width': 15, 'headwidth': 30},
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                       edgecolor="forestgreen", lw=1,))

plt.plot(temps_ana,sol_ana,"-",label='analytique',
         color="firebrick")
plt.plot(temps_exp,sol_exp,".",fillstyle='none',label='euler_exp',
         color="green")
plt.plot(temps_imp,sol_imp,".",fillstyle='none',label='euler_imp',
         color="blue")
plt.plot(temps_trap,sol_trap,".",fillstyle='none',label='trapezes',
         color="black")
plt.plot(temps_rk,sol_rk,".",fillstyle='none',label='runge_kutta',
         color="deeppink")
plt.legend(bbox_to_anchor=(0.73, 0.7), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('position',fontsize=16)
plt.xlabel('temps (secondes)',fontsize=16)
#plt.xscale('log')

plt.savefig("homogene_ordre1")
plt.show()