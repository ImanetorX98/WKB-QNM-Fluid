"""Leading scalar angular action, independent of angular eigenvalue solvers."""
import mpmath as mp


def action(A,c,mu):
    c=mp.mpf(c); mu=abs(mp.mpf(mu))
    if c==0:
        return mp.pi*(mp.sqrt(A)-mu),mp.pi/(2*mp.sqrt(A)),mp.mpf(0)
    disc=mp.sqrt((A+c*c)**2-4*c*c*mu*mu)
    tp=(c*c-A+disc)/(2*c*c)
    tm=(c*c-A-disc)/(2*c*c)
    if not (0<tp<1 and tm<0):
        raise ValueError('Outside the single equator-crossing allowed interval')
    def xx(u):return tp*mp.sin(u)**2
    I=mp.quad(lambda u:2*c*tp*mp.cos(u)**2*mp.sqrt(xx(u)-tm)/(1-xx(u)),[0,mp.pi/2])
    IA=mp.quad(lambda u:1/(c*mp.sqrt(xx(u)-tm)),[0,mp.pi/2])
    Ic=mp.quad(lambda u:2*xx(u)/mp.sqrt(xx(u)-tm),[0,mp.pi/2])
    return I,IA,Ic


def solve(c,mu):
    A=mp.findroot(lambda a:action(a,c,mu)[0]-mp.pi*(1-abs(mu)),(.85,1.))
    I,IA,Ic=action(A,c,mu)
    return A,-Ic/IA,I-mp.pi*(1-abs(mu))


if __name__=='__main__':
    for precision in (35,55):
        mp.mp.dps=precision
        mu=mp.mpf(2)/3;c=mp.mpf(3)/5
        A,derivative,res=solve(c,mu)
        step=mp.mpf('0.00001')
        fd=(solve(c+step,mu)[0]-solve(c-step,mu)[0])/(2*step)
        print('dps=',precision,'A0=',mp.nstr(A,27),'A0_c=',mp.nstr(derivative,23),
              'action_residual=',mp.nstr(abs(res),4),'FD_error=',mp.nstr(abs(fd-derivative),5))
