_author_ = "Nicolas Coruzzi"
_filename_ = "Euler_ordre1_homogene"
_creationdate_ = "01/12/19"

#x'+cx=0
#solution de la forme f(t)=Cexp(-ct)
#si pour tz on a f(tz)=xz, alors C=xz*exp(c*tz)

#yn+1= (1-c*h)yn

def euler_exp_ordre1_homogene(c,xz,t_tot,nb):
    delta_t = t_tot / nb
    x = xz
    t = 0
    X = [x]
    tps = [t]
    for n in range(nb):
        x = x + (delta_t) * (-c*x)
        t = t + delta_t
        X.append(x)
        tps.append(t)
    return (X,tps)

def euler_imp_ordre1_homogene(c,xz,t_tot,nb):
    delta_t = t_tot / nb
    x = xz
    t = 0
    X = [x]
    tps = [t]
    for n in range(nb):
        x = x /(1+delta_t*c)
        t = t + delta_t
        X.append(x)
        tps.append(t)
    return (X,tps)