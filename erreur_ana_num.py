_author_ = "Nicolas Coruzzi"
_filename_ = "erreur_ana_num"
_creationdate_ = "02/12/19"

def erreur(u,v,t_tot,nb):
    '''prend en entree deux listes, la premiere
    est celle des positions par analytique,
    et la seconde celle des positions par une methode num.
    renvoie l'erreur relative des deux mesures'''
    delta_t = t_tot / nb
    res=0
    for i in range (len(u)):
        res+=(u[i]-v[i])**2
    #on renvoie l'erreur et le delta_t
    return ((delta_t*res)**(1/2),delta_t)
