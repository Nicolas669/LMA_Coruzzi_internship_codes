_author_ = "Nicolas Coruzzi"
_filename_ = "fonctions_matrice"
_creationdate_ = "01/12/19"

import numpy as np

def inverse(matrice):
    '''uniquement pour une matrice deux deux'''
    a=matrice[0][0]
    b=matrice[0][1]
    c=matrice[1][0]
    d=matrice[1][1]
    if a*d-b*c==0:
        return "Division par zero impossible"
    else:
        return (1/(a*d-b*c))*np.array([[d,-b],[-c,a]])

def max_mod_val_prop(matrice):
    '''renvoie le module du max des valeurs propres de matrice'''
    val_prop = np.linalg.eig(matrice)[0]
    maxval = max(val_prop)
    return abs(maxval)