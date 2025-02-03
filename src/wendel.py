""" Check expressions of Wendel factorial moments """
from numpy import *

def M1(a):
    if array(a).size == 1:
        if a > 1:
            return 0
        return 1/sqrt(a) - 1
    return where(a>1, 0, 1/sqrt(a) - 1)

def M2(a):
    f = 1-2/sqrt(a)+ (1/pi)*((4/sqrt(a))*arcsin(sqrt(a/(1-a))) - 2*arcsin(a/(1-a)) + 2*sqrt(1-2*a)/a)    
    return f

def wendel_integral(a, x0,M):
    # Integral from x0 to 1-r of M(r/x) f_{r+x}(x) dx
    n = 1000000
    dx = (1-a - x0)/n
    x = x0 + (0.5 + arange(n))*dx
    
    integ = sqrt((1-a-x)/(a*x))/(pi * (1-x))
#    print('Integ in', integ.min(), integ.max(), integ[:10])
    return (M1(a/x)*integ).sum() * dx
def test_M2(a):
    print('integ', wendel_integral(a, 1*a, M1))#(M1(a/x) * integ).sum() *dx)
    print('Analytic', M2(a))


if __name__=='__main__':
    print('M1(x)', M1(0.29))
    test_M2(0.29)
    
