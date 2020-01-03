_author_ = "Nicolas Coruzzi"
_filename_ = "generateur_fonction_nonlineaire"
_creationdate_ = "18/12/19"

import matplotlib.pyplot as plt
from math import tan


raid=5.
Uz=1.
Vz=0.
u_debut=-1
u_fin=5
theta=1.5

alpha1=Uz
alpha2=3.

k1=raid
k2=-raid
k3=theta*k1

def FUNC(u,alpha1,alpha2,k1,k2,k3):
    """renvoie la valeur de la fonction [0] et de sa derivee [1] en u"""
    if u<=alpha1:
        return (k1*u,k1)
    if u>alpha1 and u<alpha2:
        return (tan(k2)*alpha2+k2*u,k2)
    if u>=alpha2:
        return ((k2-tan(k3)*alpha2)+k3*u,k3)

i=u_debut
F_res=[]
U_res=[]
while i<u_fin:
    U_res+=[i]
    F_res+=[FUNC(i,alpha1,alpha2,k1,k2,k3)[0]]
    i+=0.01

plt.figure(figsize=(14,8), dpi=80)

plt.plot(U_res,F_res,"-",fillstyle='none',label='F(u)',
         color="red")

plt.legend(bbox_to_anchor=(0.54, 0.7), loc='lower left',
           fontsize =14,borderaxespad=0.1)
plt.ylabel('F',fontsize=16)
plt.xlabel('u',fontsize=16)

plt.show()