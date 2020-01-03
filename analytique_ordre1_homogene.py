_author_ = "Nicolas Coruzzi"
_filename_ = "analytique_ordre1_homogene"
_creationdate_ = "01/12/19"

from math import exp

#x'+cx=0
#solution de la forme f(t)=Cexp(-ct)
#si pour tz on a f(tz)=xz, alors C=xz*exp(c*tz)

def analytique_ordre1_homogene_discrete(c,xz,t_dis):
    return xz*exp(-c*t_dis)

def analytique_ordre1_homogene_complete(c,xz,t_tot,nb):
    sol=[]
    tps=[]
    delta_t=t_tot/nb
    t=0
    tps+=[t]
    sol+=[xz]
    while t<=t_tot:
        t += delta_t
        tps+=[t]
        sol += [analytique_ordre1_homogene_discrete(c, xz, t)]
    return (sol,tps)