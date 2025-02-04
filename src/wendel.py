""" Check expressions of Wendel factorial moments """
from numpy import *

def Mn_raw(n, a):
    """
    Wendel factorial moments M*_n(a). This (raw) version doesnt zero values above 1/n
    """
    if n==0:
        return 1
    if n==1:
        return 1/sqrt(a) - 1
    if n==2:
        return 1-2/sqrt(a)+ (1/pi)*((4/sqrt(a))*arcsin(sqrt(a/(1-a))) - 2*arcsin(a/(1-a)) + 2*sqrt(1-2*a)/a)
    if n==3:
        """
        M1 = 3{F4-F3) + 2{F3-F2} + (1C1){F2 - F1}
        M2 = 3{F4-F3} + 1{F3-F2}
        M3 = 1{F4-F3}
        
        Prob all <a = 1 - {F2-F1} - {F3-F2} - {F4-F3]
        = 1 - (M1 - M2 + M3)
        => M3 = 1 - M1 + M2 - B
        """
        return 3/sqrt(a) - 1 + ((3*a+1)/(a*sqrt(a)) - 6*sqrt(1-2*a)/a + 6*arcsin(a/(1-a))  -12*arcsin(sqrt(a/(1-a))) / sqrt(a))/pi
    else:
        raise Exception('Not implemented')
    
def M(n,a):
    if array(a).size == 1:
        if a * n > 1:
            return 0
        return Mn_raw(n,a)
    a = array(a)
    res = zeros_like(a)
    idx = flatnonzero(n*a <= 1)
    res[idx] = Mn_raw(n, a[idx])
    return res


def wendel_integral(a, x0,M):
    # Integral from x0 to 1-r of M(r/x) f_{r+x}(x) dx
    n = 1000000
    dx = (1-a - x0)/n
    x = x0 + (0.5 + arange(n))*dx
    
    integ = sqrt((1-a-x)/(a*x))/(pi * (1-x))
#    print('Integ in', integ.min(), integ.max(), integ[:10])
    return (M(a/x)*integ).sum() * dx

def test_M(n,a):
    print('Test M%d(%f)'%(n,a))
    print('Num int.', wendel_integral(a, 0, lambda x : M(n-1,x)))
    print('Analytic', M(n,a))
    
if __name__=='__main__':
    test_M(1,0.29)
    test_M(2,0.29)
    
    test_M(3,0.26)
    test_M(3,0.16)
