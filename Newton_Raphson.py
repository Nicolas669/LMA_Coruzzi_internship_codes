_author_ = "Nicolas Coruzzi"
_filename_ = "Newton_Raphson"
_creationdate_ = "18/12/19"

import numpy as np

def fonction_cube(x):
    """renvoie la fonction cube et sa dérivee"""
    fonc=x**3 - 2*x-5
    der=3*x**2 -2
    return (fonc,der)

def newton_raphson_scalaire(precision,nb_iteration,point_depart,fonction):
    i=0
    X=point_depart
    while (abs((fonction(X)[0])/(fonction(X)[1]))>=precision) and (i<=nb_iteration):
        X=X-((fonction(X)[0])/(fonction(X)[1]))
        i+=1
    return X

def newton_raphson_vectoriel(precision,nb_iteration,point_depart,G_method,U_n,FONC):
    i=0
    X=point_depart #X est un vecteur 2x1 ici
    while (abs(np.dot((G_method(X,U_n,FONC)[0]).T,np.linalg.inv(G_method(X,U_n,FONC)[1]))).all()>=precision) and (i<=nb_iteration):
        X=(X-(np.dot((G_method(X,U_n,FONC)[0]).T,np.linalg.inv(G_method(X,U_n,FONC)[1]))))
        i+=1
    return X

#precision=10**(-10)
#nb_iteration=50
#point_depart=2

#res=newton_raphson_scalaire(precision,nb_iteration,point_depart,fonction_cube)
#print(res)
#print()
#print(fonction_cube(res)[0])