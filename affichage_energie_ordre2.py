_author_ = "Nicolas Coruzzi"
_filename_ = "affichage_energie_ordre2"
_creationdate_ = "03/12/19"

#Si on considere un coeffcient d'amortissement nul, alors l'energie de notre systeme se conserve au cours du temps

#on a  X'' + 2*omega*beta*X' + omega**2 * X = 0
#si amortissement nul, alors beta nul, et donc on a X'' + omega**2 * X = 0

#en multipliant par X', on a X'*X'' + X'*X *omega**2 + 0
#et donc on a la derivee par rapport au temps de (1/2)X'**2 + (omega**2 / 2)X**2 egale a 0
#et donc (1/2)X'**2 + (omega**2 / 2)X**2 constant --> energie se conserve


from analytique_ordre2_homogene import*
from Euler_ordre2_homogene import*
from trapezes_ordre2_homogene import*
from Runge_Kutta_4_ordre2 import*
from Newmark_ordre2_homogene import*
import matplotlib.pyplot as plt

raid=5
long=1
m=1
#ATTENTION IL FAUT UN 'ca' NUL POUR CE PROGRAMME !!!
ca=0
Xz=1
vz=0

omega = (raid / m) ** (1 / 2)

periode=(2*np.pi)/((raid / m) ** (1 / 2))

nb_periode_t_tot=5
t_tot=nb_periode_t_tot*periode
nb=50*nb_periode_t_tot

beta = ca / (m * 2 * omega)

M = np.array([[0, 1], [-(omega ** 2), -2 * omega * beta]])


(sol_exp,temps_exp)=euler_exp_ordre2_homogene(M,Xz,vz,t_tot,nb)
(sol_imp,temps_imp)=euler_imp_ordre2_homogene(M,Xz,vz,t_tot,nb)
(sol_trap,temps_trap)=trapezes_ordre2_homogene(M,Xz,vz,t_tot,nb)
(sol_rk,temps_rk)=runge_kutta_ordre2_homogene(M,Xz,vz,t_tot,nb)
(sol_nmark,temps_nmark,vitesse_nmark)=newmark_ordre2_homogene(raid,m,ca,Xz,vz,t_tot,nb)

energie_ana=[(1/2)*vz**2+(omega**2 / 2)*Xz**2]*len(temps_exp)

energie_exp=[]
energie_imp=[]
energie_trap=[]
energie_rk=[]
energie_nmark=[]

for j in range(len(temps_exp)):
    energie_exp+=[(1/2)*sol_exp[1][j]**2+(omega**2 / 2)*sol_exp[0][j]**2]
    energie_imp+=[(1/2)*sol_imp[1][j]**2+(omega**2 / 2)*sol_imp[0][j]**2]
    energie_trap+=[(1/2)*sol_trap[1][j]**2+(omega**2 / 2)*sol_trap[0][j]**2]
    energie_rk+=[(1/2)*sol_rk[1][j]**2+(omega**2 / 2)*sol_rk[0][j]**2]
    energie_nmark+=[(1/2)*vitesse_nmark[j]**2+(omega**2 / 2)*sol_nmark[j]**2]

#############Proprietes de la fenetre
plt.figure(figsize=(12,8), dpi=80)
plt.plot(temps_exp,energie_ana,"-",label='analytique',
         color="firebrick")
plt.plot(temps_exp,energie_exp,".",fillstyle='none',label='euler_exp',
         color="green")

plt.annotate(' coefficient amortissement : '+ str(ca)+
             ' \n position initiale : '+ str(Xz)+
             ' \n vitesse initiale : '+ str(vz),
             fontsize=20, xy = (0.2,0.7),
             xycoords='figure fraction', xytext = (0.2,0.7),
             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
                           'width': 15, 'headwidth': 30},
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                       edgecolor="forestgreen", lw=1,))

plt.legend(bbox_to_anchor=(0.73, 0.2), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('Energie du systeme',fontsize=16)
plt.xlabel('temps (secondes)',fontsize=16)
#plt.yscale('log')
plt.savefig("conservation_energie_ordre2_exp")
###########################################################
plt.figure(figsize=(12,8), dpi=80)
plt.plot(temps_exp,energie_ana,"-",label='analytique',
         color="firebrick")
plt.plot(temps_imp,energie_imp,".",fillstyle='none',label='euler_imp',
         color="blue")

plt.annotate(' coefficient amortissement : '+ str(ca)+
             ' \n position initiale : '+ str(Xz)+
             ' \n vitesse initiale : '+ str(vz),
             fontsize=20, xy = (0.3,0.7),
             xycoords='figure fraction', xytext = (0.3,0.7),
             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
                           'width': 15, 'headwidth': 30},
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                       edgecolor="forestgreen", lw=1,))

plt.legend(bbox_to_anchor=(0.73, 0.2), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('Energie du systeme',fontsize=16)
plt.xlabel('temps (secondes)',fontsize=16)
#plt.yscale('log')
plt.savefig("conservation_energie_ordre2_imp")
###########################################################
plt.figure(figsize=(12,8), dpi=80)
plt.plot(temps_exp,energie_ana,"-",label='analytique',
         color="firebrick")
plt.plot(temps_trap,energie_trap,".",fillstyle='none',label='trapezes',
         color="black")

plt.annotate(' coefficient amortissement : '+ str(ca)+
             ' \n position initiale : '+ str(Xz)+
             ' \n vitesse initiale : '+ str(vz),
             fontsize=20, xy = (0.2,0.2),
             xycoords='figure fraction', xytext = (0.2,0.2),
             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
                           'width': 15, 'headwidth': 30},
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                       edgecolor="forestgreen", lw=1,))

plt.legend(bbox_to_anchor=(0.73, 0.2), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('Energie du systeme',fontsize=16)
plt.xlabel('temps (secondes)',fontsize=16)
#plt.yscale('log')
plt.savefig("conservation_energie_ordre2_trap")
###########################################################
plt.figure(figsize=(12,8), dpi=80)
plt.plot(temps_exp,energie_ana,"-",label='analytique',
         color="firebrick")
plt.plot(temps_rk,energie_rk,".",fillstyle='none',label='runge_kutta',
         color="deeppink")

plt.annotate(' coefficient amortissement : '+ str(ca)+
             ' \n position initiale : '+ str(Xz)+
             ' \n vitesse initiale : '+ str(vz),
             fontsize=20, xy = (0.2,0.2),
             xycoords='figure fraction', xytext = (0.2,0.2),
             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
                           'width': 15, 'headwidth': 30},
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                       edgecolor="forestgreen", lw=1,))

plt.legend(bbox_to_anchor=(0.73, 0.2), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('Energie du systeme',fontsize=16)
plt.xlabel('temps (secondes)',fontsize=16)
#plt.yscale('log')
plt.savefig("conservation_energie_ordre2_rk")
###########################################################
plt.figure(figsize=(12,8), dpi=80)
plt.plot(temps_exp,energie_ana,"-",label='analytique',
         color="firebrick")
plt.plot(temps_nmark,energie_nmark,".",fillstyle='none',label='newmark',
         color="goldenrod")

plt.annotate(' coefficient amortissement : '+ str(ca)+
             ' \n position initiale : '+ str(Xz)+
             ' \n vitesse initiale : '+ str(vz),
             fontsize=20, xy = (0.2,0.2),
             xycoords='figure fraction', xytext = (0.2,0.2),
             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
                           'width': 15, 'headwidth': 30},
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                       edgecolor="forestgreen", lw=1,))

plt.legend(bbox_to_anchor=(0.73, 0.2), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('Energie du systeme',fontsize=16)
plt.xlabel('temps (secondes)',fontsize=16)
#plt.yscale('log')
plt.savefig("conservation_energie_ordre2_nmark")
###########################################################
plt.figure(figsize=(12,8), dpi=80)
plt.plot(temps_exp,energie_ana,"-",label='analytic',
         color="firebrick")
#plt.plot(temps_exp,energie_exp,".",fillstyle='none',label='euler_exp',
#         color="green")
plt.plot(temps_imp,energie_imp,"+",fillstyle='none',label='euler implicite',
         color="blue")
plt.plot(temps_trap,energie_trap,".",fillstyle='none',label='trapezoidal',
         color="black")
plt.plot(temps_rk,energie_rk,"x",fillstyle='none',label='Runge Kutta',
         color="deeppink")
plt.plot(temps_nmark,energie_nmark,"v",fillstyle='none',label='Newmark',
         color="goldenrod")

#plt.annotate(' coefficient amortissement : '+ str(ca)+
#             ' \n position initiale : '+ str(Xz)+
#             ' \n vitesse initiale : '+ str(vz),
#             fontsize=20, xy = (0.15,0.8),
#             xycoords='figure fraction', xytext = (0.15,0.8),
#             arrowprops = {'facecolor': 'white', 'edgecolor': 'white',
#                           'width': 15, 'headwidth': 30},
#             bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
#                       edgecolor="forestgreen", lw=1,))

plt.legend(bbox_to_anchor=(0.7, 0.5), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('System energy',fontsize=16)
plt.xlabel('time (seconds)',fontsize=16)
#plt.yscale('log')
plt.savefig("conservation_energie_ordre2_tous")

plt.show()