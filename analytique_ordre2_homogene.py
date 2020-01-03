_author_ = "Nicolas Coruzzi"
_filename_ = "analytique_ordre2_homogene"
_creationdate_ = "01/12/19"

from math import exp,cos,sin

#ces fonctions sont pour notre equation differentielle pour le ressort seul

def analytique_ordre2_homogene_discrete(omega,beta,Xz,vz,t_dis):

    if beta>1:
        return ((((beta-((beta**2)-1)**(1/2))*Xz+(vz/omega))/(-2*((beta**2)-1)**(1/2)))*exp((-beta-((beta**2)-1)**(1/2))*omega*t_dis) + ((-(beta+((beta**2)-1)**(1/2))*Xz-(vz/omega))/(-2*((beta**2)-1)**(1/2)))*exp((-beta+((beta**2)-1)**(1/2))*omega*t_dis))
    if beta<1:
        return ((Xz*cos(omega*(1-beta**2)**(1/2)*t_dis)+(((vz/omega)+beta*Xz)/(1-beta**2)**(1/2))*sin(omega*(1-beta**2)**(1/2)*t_dis))*exp(-omega*beta*t_dis))
    if beta==1:
        return (((vz+omega*beta*Xz)*t_dis+Xz)*exp(-omega*beta*t_dis))

def analytique_ordre2_homogene_complete(omega,beta,Xz,vz,t_tot,nb):
    sol=[]
    tps=[]
    delta_t=t_tot/nb
    sol+=[Xz]
    tn = 0
    tps += [tn]
    for n in range(nb):
        tn = tn + delta_t
        tps += [tn]
        sol += [analytique_ordre2_homogene_discrete(omega,beta, Xz, vz, tn)]
    return (sol,tps)