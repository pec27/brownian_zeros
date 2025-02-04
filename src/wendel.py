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
    f = 1-1/sqrt(a)+ (1/pi)*(-(2/sqrt(a))*arcsin((1-3*a)/(1-a)) - 2*arcsin(a/(1-a)) + 2*sqrt(1-2*a)/a)    

    return f

def wendel_integral(a, x0,M):
    # Integral from x0 to 1-r of M(r/x) f_{r+x}(x) dx
    n = 1000000
    dx = (1-a - x0)/n
    x = x0 + (0.5 + arange(n))*dx
    
    integ = sqrt((1-a-x)/(a*x))/(pi * (1-x))
#    print('Integ in', integ.min(), integ.max(), integ[:10])
    return (M(a/x)*integ).sum() * dx

def test_M1(a):
    print('Num int.', wendel_integral(a, 0, lambda x : 1))
    print('Analytic', M1(a))
def test_M2(a):
    print('Num int.', wendel_integral(a, 1*a, M1))#(M1(a/x) * integ).sum() *dx)
    print('Analytic', M2(a))

def M3(a):
    """
    M1 = 3{F4-F3) + 2{F3-F2} + (1C1){F2 - F1}
    M2 = 3{F4-F3} + 1{F3-F2}
    M3 = 1{F4-F3}
    
    Prob all <a = 1 - {F2-F1} - {F3-F2} - {F4-F3]
    = 1 - (M1 - M2 + M3)
    => M3 = 1 - M1 + M2 - B
    
    """
    # M3 = P3 + 3 P2 + 3P1 + P0, where P0 = F1 - F0
    # F3 = P3 + P2 + P1 + P0
    #    = M3 - 2M2 + 2M1
    # => M3 = F3 + 2M2 - 2M1
    # The way we have defined it, B is the probability that the largest (and hence all) zero-free interval < a.
    # In Wendel, F_N(x) is the probability that the Nth interval has length < a, hence
    # F_{N+1}(x) - F_N(x) = (1-F_N)-(1-F_{N+1}) is the probability that the Nth interval is longer than a, but the N+1th is shorter,
    # hence there are N intervals longer than a. In order for this to work for N=0, F_0(a)=0 for all a.
    # Now we know B is the probability that any interval longer than a, so B = 1-F_1. 
    B =  (3*a+1)/(a*sqrt(a)) - 6*sqrt(1-2*a)/a + 6*arcsin(a/(1-a)) + 6*arcsin((1-3*a)/(1-a)) / sqrt(a)
    B =  (3*a+1)/(a*sqrt(a)) - 6*sqrt(1-2*a)/a + 6*arctan(a/sqrt(1-2*a)) + 6*arctan((1-3*a)/sqrt(-8*a*a+4*a)) / sqrt(a)
    B =  (3*a+1)/(a*sqrt(a)) - 6*sqrt(1-2*a)/a + 6*arcsin(a/(1-a))  -12*arcsin(sqrt(a/(1-a))) / sqrt(a)

    B = 3/sqrt(a) + (1/pi)*B - 1
    
     
    R = B
    return R
    
def test_M3(a):
    print('Num. int', wendel_integral(a, 2*a, M2))

    print('Analytic', M3(a))


if __name__=='__main__':
    #print('M1(x)', M1(0.29))
    test_M2(0.29)
    
    test_M3(0.26)
    test_M3(0.16)
    #test_M1(0.3)
