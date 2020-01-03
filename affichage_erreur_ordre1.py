_author_ = "Nicolas Coruzzi"
_filename_ = "affichage_erreur_ordre1"
_creationdate_ = "02/12/19"

from analytique_ordre1_homogene import*
from Euler_ordre1_homogene import*
from trapezes_ordre1_homogene import*
from Runge_Kutta_4_ordre1 import*
import matplotlib.pyplot as plt
import numpy as np
from erreur_ana_num import*
from math import atan

c=5
xz=10
t_tot=5
nb=[50,100,200,300,400,500,600]

res_exp=[]
res_imp=[]
res_trap=[]
res_rk=[]
delta_t_exp=[]
delta_t_imp=[]
delta_t_trap=[]
delta_t_rk=[]
for nombre in nb:
    (sol_ana,temps_ana)=analytique_ordre1_homogene_complete(c,xz,t_tot,nombre)
    (sol_exp,temps_exp)=euler_exp_ordre1_homogene(c,xz,t_tot,nombre)
    (sol_imp,temps_imp)=euler_imp_ordre1_homogene(c,xz,t_tot,nombre)
    (sol_trap,temps_trap)=trapezes_ordre1_homogene(c,xz,t_tot,nombre)
    (sol_rk, temps_rk) = runge_kutta_ordre1_homogene(c, xz, t_tot, nombre)

    (err_exp,delta_exp)=erreur(sol_ana,sol_exp,t_tot,nombre)
    (err_imp,delta_imp)=erreur(sol_ana,sol_imp,t_tot,nombre)
    (err_trap,delta_trap)=erreur(sol_ana,sol_trap,t_tot,nombre)
    (err_rk, delta_rk) = erreur(sol_ana, sol_rk, t_tot, nombre)

    res_exp+=[err_exp]
    res_imp+=[err_imp]
    res_trap+=[err_trap]
    res_rk+=[err_rk]
    delta_t_exp+=[delta_exp]
    delta_t_imp+=[delta_imp]
    delta_t_trap+=[delta_trap]
    delta_t_rk+=[delta_rk]

def calcul_pente_gauche(s,t):
    '''prend en entree 2 listes et renvoie la pente entre les deux plus petits points de ces listes'''

    for i in range(len(s)):
        if s[i]==min(s) and t[i]==min(t):
            a=t[i]
            fa=s[i]
            del s[i]
            del t[i]
            for j in range(len(s)):
                if s[j]==min(s) and t[j]==min(t):
                    b=t[j]
                    fb=s[j]
                    return atan((fb-fa)/(b-a))

#############Proprietes de la fenetre
pente_trap=calcul_pente_gauche(res_trap,delta_t_trap)
#res_pente_trap=[]
#for elt_trap in delta_t_trap:
#    res_pente_trap+=[pente_trap*elt_trap]

plt.figure(figsize=(12,8), dpi=80)
plt.plot(delta_t_trap,res_trap,"o-",label='erreur_trapezes',
         color="black")

plt.annotate(' pente a gauche : '+ str(pente_trap),
             fontsize=20, xy = (0.2,0.92),
             xycoords='figure fraction', xytext = (0.2,0.92),
             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
                           'width': 15, 'headwidth': 30},
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                       edgecolor="forestgreen", lw=1,))

#plt.plot(delta_t_trap,res_pente_trap,"--",label='pente:'+str(pente_trap),
#         color="grey")

plt.legend(bbox_to_anchor=(0.73, 0.2), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('erreur',fontsize=16)
plt.xlabel('delta_t',fontsize=16)
plt.xscale('log')
plt.yscale('log')

plt.savefig("erreur_ordre1_trapeze")
###########################################################
pente_imp=calcul_pente_gauche(res_imp,delta_t_imp)
#res_pente_imp=[]
#for elt_imp in delta_t_imp:
#    res_pente_imp+=[pente_imp*elt_imp]

plt.figure(figsize=(12,8), dpi=80)
plt.plot(delta_t_imp,res_imp,"o-",label='erreur_euler_imp',
         color="blue")

plt.annotate(' pente a gauche : '+ str(pente_imp),
             fontsize=20, xy = (0.2,0.92),
             xycoords='figure fraction', xytext = (0.2,0.92),
             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
                           'width': 15, 'headwidth': 30},
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                       edgecolor="forestgreen", lw=1,))

#plt.plot(delta_t_imp,res_pente_imp,"--",label='pente:'+str(pente_imp),
#         color="grey")
plt.legend(bbox_to_anchor=(0.73, 0.2), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('erreur',fontsize=16)
plt.xlabel('delta_t',fontsize=16)
plt.xscale('log')
plt.yscale('log')

plt.savefig("erreur_ordre1_implicite")
###########################################################
pente_exp=calcul_pente_gauche(res_exp,delta_t_exp)
#res_pente_exp=[]
#for elt_exp in delta_t_exp:
#    res_pente_exp+=[pente_exp*elt_exp]

plt.figure(figsize=(12,8), dpi=80)
plt.plot(delta_t_exp,res_exp,"o-",label='erreur_euler_exp',
        color="green")

plt.annotate(' pente a gauche : '+ str(pente_exp),
             fontsize=20, xy = (0.2,0.92),
             xycoords='figure fraction', xytext = (0.2,0.92),
             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
                           'width': 15, 'headwidth': 30},
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                       edgecolor="forestgreen", lw=1,))

#plt.plot(delta_t_exp,res_pente_exp,"--",label='pente:'+str(pente_exp),
#         color="grey")
plt.legend(bbox_to_anchor=(0.73, 0.2), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('erreur',fontsize=16)
plt.xlabel('delta_t',fontsize=16)
plt.xscale('log')
plt.yscale('log')

plt.savefig("erreur_ordre1_explicite")
###########################################################
pente_rk=calcul_pente_gauche(res_rk,delta_t_rk)
#res_pente_rk=[]
#for elt_rk in delta_t_rk:
#    res_pente_rk+=[pente_rk*elt_rk]

plt.figure(figsize=(12,8), dpi=80)
plt.plot(delta_t_rk,res_rk,"o-",label='erreur_runge_kutta',
        color="deeppink")

plt.annotate(' pente a gauche : '+ str(pente_rk),
             fontsize=20, xy = (0.2,0.92),
             xycoords='figure fraction', xytext = (0.2,0.92),
             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
                           'width': 15, 'headwidth': 30},
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                       edgecolor="forestgreen", lw=1,))

#plt.plot(delta_t_rk,res_pente_rk,"--",label='pente:'+str(pente_rk),
#         color="grey")
plt.legend(bbox_to_anchor=(0.73, 0.2), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('erreur',fontsize=16)
plt.xlabel('delta_t',fontsize=16)
plt.xscale('log')
plt.yscale('log')

plt.savefig("erreur_ordre1_runge_kutta")


plt.show()