_author_ = "Nicolas Coruzzi"
_filename_ = "Newmark_ordre2_homogene"
_creationdate_ = "03/12/19"

def newmark_ordre2_homogene(raid,m,ca,Xz,vz,t_tot,nb):
    #equation differentielle X'' + 2*omega*beta*X' + omega**2*X = 0

    #dans ce programme, D est la position, V la vitesse et A l'acceleration

    tps = []
    delta_t = t_tot / nb

    M=1

    #C=2*omega*beta (mais pas le beta de newmark!)
    C=2*((raid / m) ** (1 / 2))*(ca / (m * 2 * ((raid / m) ** (1 / 2))))

    #K=omega**2
    K=(raid / m)

    gamma=1/2
    beta=1/4

    D=[Xz]
    V=[vz]
    A=[-C*V[0]-K*D[0]]

    S = M + gamma*delta_t*C + beta*(delta_t**2)*K

    tn=0
    tps+=[tn]
    for n in range(1,nb+1):
        tn = tn + delta_t
        tps+=[tn]

        Dpr= D[n-1]+delta_t*V[n-1]+(delta_t**2)*((1/2)-beta)*A[n-1]
        Vpr= V[n-1]+delta_t*(1-gamma)*A[n-1]

        theta= -C*Vpr -K*Dpr

        A+=[theta/S]
        D+=[Dpr+beta*(delta_t**2)*A[n]]
        V+=[Vpr+gamma*delta_t*A[n]]
    return (D,tps,V)