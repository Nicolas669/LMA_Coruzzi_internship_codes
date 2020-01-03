_author_ = "Nicolas Coruzzi"
_filename_ = "trapezes_ordre1_homogene"
_creationdate_ = "01/12/19"

#x'+cx=0
#solution de la forme f(t)=Cexp(-ct)
#si pour tz on a f(tz)=xz, alors C=xz*exp(c*tz)

#yn+1= yn +h/2( f(xn,yn) + f(xn+1, yn+1))

def trapezes_ordre1_homogene(c,xz,t_tot,nb):
    delta_t = t_tot / nb
    x = xz
    t = 0
    X = [x]
    tps = [t]
    for n in range(nb):
        x = x *((1-(c*delta_t)/2)/(1+(c*delta_t)/2))
        t = t + delta_t
        X.append(x)
        tps.append(t)
    return (X,tps)