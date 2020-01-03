_author_ = "Nicolas Coruzzi"
_filename_ = "affichage_Nmasses_nonlineaire"
_creationdate_ = "20/12/19"

import matplotlib.pyplot as plt

raid=5.
m=1.
ca=0.
Utz=1.
Vtz=0.
#Utz et Vtz sont ici la position et vitesse de la masse n+1 a l'instant initial (on suppose les autres au repos)

N=2

nb=500
t_tot=10

#pour newton raphson
precision=10**(-10)
nb_iteration=10
####################
import numpy as np

Y = np.zeros((2 * N, nb + 1), float)
Y[(2 * N - 2)][0] = Utz
Y[(2 * N - 1)][0] = Vtz

def fonction(x,k):
    """renvoie la fonction et sa dérivee"""
    fon=k*(x+x**3)
    der=(1+3*(x**2))*k
    return (fon,der)


#Attention la FONCTION n'est pas bonne encore ici!!!
def FONCTION(elt,alph,m,fonc,k):

    F=np.zeros((2*elt.shape[0],1),float)

    F[0]=elt[1]
    F[1]=(-alph / m) * elt[1]-((fonc(elt[0],k)[0]) / m)-((fonc((elt[0]-elt[2]),k)[0]) / m)
    for i in range(2,elt.shape[0]-3):
        F[i]=elt[i+1]
        F[i+1]=(-alph / m) * elt[i+1]-((fonc((elt[i]-elt[i+2]),k)[0]) / m)-((fonc((elt[i]-elt[i-2]),k)[0]) / m)
    F[2*elt.shape[0]-1]=elt[2*elt.shape[0]-1]
    F[2*elt.shape[0]]=(-alph / m) * elt[2*elt.shape[0]-1]-((fonc((elt[2*elt.shape[0]-1]-elt[2*elt.shape[0]-2]),k)[0]) / m)-((fonc((elt[i]-elt[i-2]),k)[0]) / m)

    F1 = elt[1]
    F2 = ((-alph / m) * elt[1]) - ((fonc(elt[0],k)[0]) / m)
    F = np.array([F1, F2])
    JacF = np.array([[0., 1.], [-((fonc(elt[0],k)[1]) / m), (-alph / m)]])
    return (F, JacF)

print(Y[:,0].shape[0])
print()
Y[:,0][1]=5
print(Y[:,0])
