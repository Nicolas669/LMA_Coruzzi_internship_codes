_author_ = "Nicolas Coruzzi"
_filename_ = "affichage_erreur_ordre2"
_creationdate_ = "02/12/19"

from analytique_ordre2_homogene import*
from Euler_ordre2_homogene import*
from trapezes_ordre2_homogene import*
from Runge_Kutta_4_ordre2 import*
from Newmark_ordre2_homogene import*
import matplotlib.pyplot as plt
import numpy as np
from erreur_ana_num import*
from math import atan,log
from copy import deepcopy

def calcul_toutes_pentes(s,t):
    res=[]
    for i in range(0,len(s)-1,1):
        b=t[i]
        fb=s[i]
        a=t[i+1]
        fa=s[i+1]
        res+=[atan((fb-fa)/(b-a))]
    return res

raid=1
long=1
m=1
ca=0
Xz=1.
vz=0.

omega = (raid / m) ** (1 / 2)
beta = ca / (m * 2 * omega)

M = np.array([[0, 1], [-(omega ** 2), -2 * omega * beta]])

periode=(2*np.pi)/((raid / m) ** (1 / 2))

t_tot=12.56

nb=[32,64,128,256,512,1024]


res_exp=[]
res_imp=[]
res_trap=[]
res_rk=[]
res_nmark=[]
delta_t_exp=[]
delta_t_imp=[]
delta_t_trap=[]
delta_t_rk=[]
delta_t_nmark=[]

res_exp_log=[]
res_imp_log=[]
res_trap_log=[]
res_rk_log=[]
res_nmark_log=[]
delta_t_exp_log=[]
delta_t_imp_log=[]
delta_t_trap_log=[]
delta_t_rk_log=[]
delta_t_nmark_log=[]

for nombre in nb:
    (sol_ana, temps_ana) = analytique_ordre2_homogene_complete(omega,beta, Xz, vz, t_tot, nombre)
    (sol_exp, temps_exp) = euler_exp_ordre2_homogene(M, Xz, vz, t_tot, nombre)
    (sol_imp, temps_imp) = euler_imp_ordre2_homogene(M, Xz, vz, t_tot, nombre)
    (sol_trap, temps_trap) = trapezes_ordre2_homogene(M, Xz, vz, t_tot, nombre)
    (sol_rk, temps_rk) = runge_kutta_ordre2_homogene(M, Xz, vz, t_tot, nombre)
    (sol_nmark, temps_nmark,vitesse_nmark) = newmark_ordre2_homogene(raid, m, ca, Xz, vz, t_tot, nombre)

    (err_exp,delta_exp)=erreur(sol_ana,sol_exp[0],t_tot,nombre)
    (err_imp,delta_imp)=erreur(sol_ana,sol_imp[0],t_tot,nombre)
    (err_trap,delta_trap)=erreur(sol_ana,sol_trap[0],t_tot,nombre)
    (err_rk, delta_rk) = erreur(sol_ana, sol_rk[0], t_tot, nombre)
    (err_nmark, delta_nmark) = erreur(sol_ana, sol_nmark, t_tot, nombre)

    res_exp_log+=[log(err_exp)]
    res_imp_log+=[log(err_imp)]
    res_trap_log+=[log(err_trap)]
    res_rk_log += [log(err_rk)]
    res_nmark_log += [log(err_nmark)]
    delta_t_exp_log+=[log(delta_exp)]
    delta_t_imp_log+=[log(delta_imp)]
    delta_t_trap_log+=[log(delta_trap)]
    delta_t_rk_log+=[log(delta_rk)]
    delta_t_nmark_log += [log(delta_nmark)]

    res_exp += [err_exp]
    res_imp += [err_imp]
    res_trap += [err_trap]
    res_rk += [err_rk]
    res_nmark += [err_nmark]
    delta_t_exp += [delta_exp]
    delta_t_imp += [delta_imp]
    delta_t_trap += [delta_trap]
    delta_t_rk += [delta_rk]
    delta_t_nmark += [delta_nmark]

for nombre in range(len(nb)):

    ordres_exp = calcul_toutes_pentes(res_exp_log, delta_t_exp_log)
    ordres_imp=calcul_toutes_pentes(res_imp_log, delta_t_imp_log)
    ordres_trap = calcul_toutes_pentes(res_trap_log, delta_t_trap_log)
    ordres_rk=calcul_toutes_pentes(res_rk_log, delta_t_rk_log)
    ordres_nmark=calcul_toutes_pentes(res_nmark_log, delta_t_nmark_log)

    if nombre != 0:
        print(nb[nombre])
        print()
        print("explicite : erreur " + str(res_exp[nombre]) + " || ordre : "+ str(ordres_exp[nombre-1]))
        print("implicite : erreur " + str(res_imp[nombre]) + " || ordre : "+ str(ordres_imp[nombre-1]))
        print("trapezes : erreur " + str(res_trap[nombre]) + " || ordre : "+ str(ordres_trap[nombre-1]))
        print("runge_kutta : erreur " + str(res_rk[nombre]) + " || ordre : "+ str(ordres_rk[nombre-1]))
        print("newmark : erreur " + str(res_nmark[nombre]) + " || ordre : "+ str(ordres_nmark[nombre-1]))
        print()

    else:
        print(nb[nombre])
        print()
        print("explicite : erreur " + str(res_exp[nombre]))
        print("implicite : erreur " + str(res_imp[nombre]))
        print("trapezes : erreur " + str(res_trap[nombre]))
        print("runge_kutta : erreur " + str(res_rk[nombre]))
        print("newmark : erreur " + str(res_nmark[nombre]))
        print()

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


################################Calcul pente trap autre methode:
#rs_trap=np.array(deepcopy(res_trap))
#delt_t_trap=np.array(deepcopy(delta_t_trap))
#print(delt_t_trap)
#a = min(delt_t_trap)
#print(a)
#ind1 = np.where(delt_t_trap == a)[0][0]
#print(ind1)
#delt_t_trap = np.delete(delt_t_trap, ind1)
#fa = rs_trap[ind1]
#rs_trap = np.delete(rs_trap, ind1)
#b = min(delt_t_trap)
#ind2 = np.where(delt_t_trap == b)[0][0]
#fb = rs_trap[ind2]

#pente_trap2= atan((fb - fa) / (b - a))
#print(pente_trap2)
###############################################################

#############Proprietes de la fenetre
pente_trap=calcul_pente_gauche(res_trap_log,delta_t_trap_log)

plt.figure(figsize=(12,8), dpi=80)
plt.plot(delta_t_trap_log,res_trap_log,"o-",label='trapezoidal error',
         color="black")

plt.legend(bbox_to_anchor=(0.73, 0.2), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('log error',fontsize=16)
plt.xlabel('log delta_t',fontsize=16)

plt.savefig("erreur_ordre2_trapeze")
###########################################################
pente_imp=calcul_pente_gauche(res_imp_log,delta_t_imp_log)

plt.figure(figsize=(12,8), dpi=80)
plt.plot(delta_t_imp_log,res_imp_log,"o-",label='euler implicite error',
         color="blue")


plt.legend(bbox_to_anchor=(0.7, 0.2), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('log error',fontsize=16)
plt.xlabel('log delta_t',fontsize=16)

plt.savefig("erreur_ordre2_implicite")
###########################################################
pente_exp=calcul_pente_gauche(res_exp_log,delta_t_exp_log)

plt.figure(figsize=(12,8), dpi=80)
plt.plot(delta_t_exp_log,res_exp_log,"o-",label='erreur_euler_exp',
        color="green")


plt.legend(bbox_to_anchor=(0.73, 0.2), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('log error',fontsize=16)
plt.xlabel('log delta_t',fontsize=16)

plt.savefig("erreur_ordre2_explicite")
###########################################################
pente_rk=calcul_pente_gauche(res_rk_log,delta_t_rk_log)

plt.figure(figsize=(12,8), dpi=80)
plt.plot(delta_t_rk_log,res_rk_log,"o-",label='Runge Kutta error',
        color="deeppink")


plt.legend(bbox_to_anchor=(0.73, 0.2), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('log error',fontsize=16)
plt.xlabel('log delta_t',fontsize=16)

plt.savefig("erreur_ordre2_runge_kutta")
###########################################################
pente_nmark=calcul_pente_gauche(res_nmark_log,delta_t_nmark_log)

plt.figure(figsize=(12,8), dpi=80)
plt.plot(delta_t_nmark_log,res_nmark_log,"o-",label='Newmark error',
        color="goldenrod")


plt.legend(bbox_to_anchor=(0.73, 0.2), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('log error',fontsize=16)
plt.xlabel('log delta_t',fontsize=16)

plt.savefig("erreur_ordre2_newmark")

plt.show()