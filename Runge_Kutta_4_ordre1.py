_author_ = "Nicolas Coruzzi"
_filename_ = "Runge_Kutta_4_ordre1"
_creationdate_ = "02/12/19"

#x'+cx=0
#solution de la forme f(t)=Cexp(-ct)
#si pour tz on a f(tz)=xz, alors C=xz*exp(c*tz)

def runge_kutta_ordre1_homogene(c,xz,t_tot,nb):
    delta_t = t_tot / nb
    x = xz
    t = 0
    X = [x]
    tps = [t]
    for n in range(nb):
        x = x + (delta_t/6)*(-c*x+2*(-c*(1-(c*delta_t)/2)*x)+2*(-c*(1-(c*delta_t)/2+((c*delta_t)/2)**2)*x)-c*((1-c*delta_t+((c*delta_t)**2)/2-(((c*delta_t)**3)/4))*x))
        t = t + delta_t
        X.append(x)
        tps.append(t)
    return (X,tps)