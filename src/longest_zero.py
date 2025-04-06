""" 
Calculations for the CDF of the longest interval between zeros 
Copyright (C) 2024 Peter Creasey
"""
# TODO use sampling formula to make numerical_LR version that is more accurate
from numpy import *

def numerical_LR(r, nm=5000,np=5003):
    """ Calculate the LR-segment, where the central interval is shorter than r numerically """
    assert(r>1/3 and r<0.5)
    # Values 0-0.5 (not including ends)
    dtm = r / nm
    tm = (0.5-r) + (0.5 + arange(nm)) * dtm
    # Values 0.5-1 (not including ends)
    dtp = r / np
    tp = 0.5 + (0.5 + arange(np)) * dtp


    # Prob that left is > r
    p_left = where(tm > r, sqrt(tm/r)-1, 0)
    # Prob that right is > r
    p_right = where(1-tp > r, sqrt((1-tp)/r)-1, 0)

    # Prob either
    p_lr = 1-outer(1-p_left, 1-p_right)
    # Width of central
    c = add.outer(-tm, tp)
    not_c = where(c < r, 1, 0)

    # Prob dens of t-,t+
    rho = 1/(2*pi * sqrt(outer(tm,1-tp)*(c*c*c)))
    return (not_c * p_lr * rho).sum() * dtm * dtp

def get_tp_tm(np=2003, nm=2000, u1_min=0):
    """ t+ shape (np,) and t- shape (nm,mp) of the given distribution """
    du1 = (1.0 - u1_min) / np
    u1 = u1_min + (0.5 + arange(np)) * du1 
    tp = 1/(1+square(sin(0.5*pi*u1)))
    

    # Values t- sampled given t+
    du2 = 1.0 / nm
    u2_2 = square((0.5 + arange(nm)) * du2)
    tm = outer(u2_2, tp) / add.outer(u2_2, 2*tp-1)
    return tp,tm,du1,du2

def numerical_cdf_LR(r, nm=2000,np=2003):
    """ Like numerical_LR, but uniformly sampling from the CDF (see paper) """
    assert(r>1/3 and r<0.5)
    # Values u1 0-1 (not including ends)
    # ignore this since we also wanted to sample C^1/2_x
    
    #u1_min = arcsin(sqrt(1/(0.5 + r)-1))*2.0/pi
#    print('u1 min', u1_min)
    u1_min = 0
    tp, tm,du1,du2 = get_tp_tm(np,nm,u1_min)
    
    # Prob that left is > r
    p_left = where(tm > r, sqrt(tm/r)-1, 0)
    # Prob that right is > r
    p_right = where(1-tp > r, sqrt((1-tp)/r)-1, 0)

    # Prob either
    p_lr = 1-(1-p_left)*(1-p_right)
    # Width of central
    c = tp-tm
    not_c = where(c < r, 1, 0)
    c_x = where(logical_and(c>r, c<0.5), 1, 0)
    
    # No need for prob dens since we are sampling from CDF
    lr = (not_c * p_lr).sum() * du1*du2
    c_x = c_x.sum() * du1*du2
    return c_x, lr

def half_numerical_LR(r, np=2000000):
    # t+ goes from 0.5 to 1-r
    dtp = (0.5-r) / np
    tp = 0.5 + (0.5 + arange(np)) * dtp
    # integrand
#    f = (2/sqrt(r*(1-tp)) - 1/r)*(1/sqrt(tp-0.5) - 1/sqrt(tp-r))/pi + (1/sqrt(r) - 3/sqrt(1-tp))*(2/sqrt(2*tp-1) - 1/sqrt(r*(tp-r)))/(pi*tp)

    f = (1/(pi * tp))*(2*sqrt((1-tp)/r) - 2)*(sqrt(0.5/((1-tp)*(tp-0.5))) - sqrt((tp-r)/((1-tp)*r))) - \
        (1/(pi*tp))*(sqrt(0.5/((1-tp)*(tp-0.5))) - sqrt(r/((1-tp)*(tp-r)))) + \
        (1/(pi*sqrt(r)))*(2-sqrt((1-tp)/r))*(1/sqrt((1-tp)*(tp-0.5)) - 1/sqrt((1-tp)*(tp-r)))
    
    LR = f.sum() * dtp
    return LR

def integral_LR(r):
    # t+ goes from 0.5 to 1-r
    tp0, tp1 = 0.5, 1-r

    # Analytic cpt (indefinite integrals)
    lim = lambda x : (4/(pi * sqrt(r)))*(arctan(sqrt(2*x-1)) - sqrt((x-r)/r) + arctan(sqrt((x-r)/r))) +\
        (-2/(pi))*(2*arctan(sqrt((2*x-1)/(1-x))) + 2*arctan(sqrt((x-r)/(r*(1-x)))) -(2/sqrt(r))*arcsin(sqrt((x-r)/(1-r)))) + \
    (-2/(pi))*(arctan(sqrt((2*x-1)/(1-x))) - arctan(sqrt((x-r)/(r*(1-x))))) +\
        (-4/(pi*sqrt(r)))*(arcsin(sqrt(2-2*x)) + arctan(sqrt((x-r)/(1-x)))) +\
        (1/(pi*r))*(-sqrt(4*x-2)+2*sqrt(x-r)) 

    ana_integ = lim(tp1)-lim(tp0)
    print('  indefinite cpt', ana_integ)
    # Analytic cpt (indefinite integrals)
    lim = lambda x : (4/(pi * sqrt(r)))*(arccos(1/sqrt(2*x)) + arccos(sqrt(r/x)) -arccos(sqrt(2*x-1))) +\
        (-2/(pi))*(3*arccos(sqrt((1-x)/x)) + arcsin(sqrt((x-r)/(x*(1-r))))) +\
        (1/(pi*r))*(-sqrt(4*x-2)-2*sqrt(x-r)) 

    ana_integ = lim(tp1)-lim(tp0)
    print('  indefinite cpt', ana_integ)

    print('  indefinite cpt', ana_integ)
    # Analytic cpt (indefinite integrals)
    lim = lambda x : (4/(pi * sqrt(r)))*(arccos(1/sqrt(2*x)) + arccos(sqrt(r/x)) -arccos(sqrt(2*x-1))) +\
        (-2/(pi))*(3*arccos(sqrt((1-x)/x)) + arcsin(sqrt((x-r)/(x*(1-r))))) +\
        (1/(pi*r))*(-sqrt(4*x-2)-2*sqrt(x-r)) 

    lim1 =  (2/pi)*(1/sqrt(r)- 1)*(arccos((3*r-1)/(1-r)) + arccos(r/(1-r))) +\
        (1/(pi*r))*(-2*sqrt(1-2*r)) 

    ana_integ = lim1
    print('  indefinite cpt', ana_integ)
    
    return ana_integ

def integral_C_half_x(r):
    ana_integ = (2/(pi*sqrt(r)))*arcsin(r/(1-r)) + (4/pi)*arcsin(sqrt((1-2*r)/(1-r))) - sqrt(2)
    ana_integ = (2/(pi*sqrt(r)))*arcsin(r/(1-r)) + (2/pi)*arccos((3*r-1)/(1-r)) - sqrt(2)    
    return ana_integ

def plot():
    nr = 500
    r_min = 0.2 #1/3
    dr = (1.0-r_min) / nr

    r = r_min + arange(nr+1)*dr
    tot_sum = where(r<0.5, (2/pi)*((1/sqrt(r))*arccos((3*r-1)/(1-r)) - arccos(r/(1-r))) +\
                    (1/(pi*r))*(-2*sqrt(1-2*r)) + 1/sqrt(r)  - 1, sqrt(1/r) - 1)    
    dens = -diff(tot_sum)/dr
    # Mid points
    r = r_min + (0.5 + arange(nr))*dr
    import pylab as pl
    pl.ylim(0,2)
    pl.plot(r,dens)
    pl.show()

def test_F2_plot():
    print("F2 (TODO currently complementary of) the CDF between 1/3 and 1/2")
    plot() 
    
    r = 0.35 # Expression is for values between 1/3 and 1/2
    
    print('r', r)
    print('Numerical', numerical_LR(r))
    num_cx, num_lr = numerical_cdf_LR(r)
    print('Numerical from CDF', num_lr)
    
    print('Half-numerical',  half_numerical_LR(r))
    LR = integral_LR(r)
    print('Integral',  LR)

    print('Numerical C^1/2_x', num_cx)
    C_x = integral_C_half_x(r)
    print('Analytic C^1/2_x', C_x)

    # Needs to include C^1_1/2 = sqrt(2) -1 (I think)
    total = C_x + LR
    tot_sum = (2/pi)*((1/sqrt(r))*arccos((3*r-1)/(1-r)) - arccos(r/(1-r))) +\
        (1/(pi*r))*(-2*sqrt(1-2*r)) + 1/sqrt(r)  - sqrt(2)    
    print('Total', total, tot_sum)

def numerical_CDF_F3(r, np,nm):
    tp, tm,du1,du2 = get_tp_tm(np,nm,0)
    # Prob that left is < r
    p_left = where(tm > r, 2-sqrt(tm/r), 1)
    # Prob that right is < r
    p_right = where(1-tp > r, 2-sqrt((1-tp)/r), 1)

    # Prob both
    p_lr = p_left*p_right
    # Width of central
    c = tp-tm
    c = where(c < r, 1, 0)

    # extra calc for right triangle F0(r/t-)F1(r/(1-tp))
    R = where(tm<r, p_right,0)*c
    print('Right triangle', R.sum()*du1*du2)
    # extra calc for lower trapezium
    TP = where(logical_and(logical_and(tm>r, tp<2*r), 1-tp>tm), p_lr, 0).sum()*du1*du2
    print('Trapezium', TP)
    # No need for prob dens since we are sampling from CDF
    F3 = (c * p_lr).sum() * du1*du2
    return F3

def numerical_CDF_F2(r, np,nm):
    tp, tm,du1,du2 = get_tp_tm(np,nm,0)
    # Prob that left is < r
    p_left = where(tm > r, 2-sqrt(tm/r), 1)
    # Prob that right is < r
    p_right = where(1-tp > r, 2-sqrt((1-tp)/r), 1)

    # Prob both
    p_lr = p_left*p_right
    # Width of central
    c = tp-tm
    c = where(c < r, 1, 0)

    # No need for prob dens since we are sampling from CDF
    F2 = (c * p_lr).sum() * du1*du2
    return F2

def integ(f, x_min,x_max, n=2000000):
    dx = (x_max-x_min)/n
    x = x_min + (0.5 + arange(n))*dx
    return f(x).sum() * dx

def F1(r):
    if array(r).size == 1:
        assert(r>= 0.5)
        assert(r <= 1)
    else:
        assert(min(r)>=0.5)
        assert(max(r)<= 1)
    return 2 - sqrt(1/r)


def F3(r):

    # rho(x,y) = 1/pi ((1-2x)(1-2y)(x+y)^3)^-1/2
    # Integral of rho w.r.t. x
    def R1(x):
        indef = lambda y : (-2/(pi*(2*y+1))) * sqrt((1-2*x)/((1-2*y)*(x+y)))
        return indef
    def R2(x):
        indef = lambda y : (-2/pi) / sqrt((1-2*y)*(x+y))
        return indef
        
    # 
    v = lambda y : (R1(r-y)(y) - R1(0.5-r)(y)) * F1(r/(0.5-y))
    tri1 = integ(v, 0,2*r-0.5)
    print(tri1)
    S = lambda y : F1(r/(0.5-y)) * (2*(R1(0.5-r)(y) -R1(y)(y)) - (1/sqrt(2*r))*(R2(0.5-r)(y) - R2(y)(y)))
    S = integ(S ,0, 2*r-0.5)
    S2 = lambda y : F1(r/(0.5-y)) * (2*(R1(r-y)(y) -R1(y)(y)) - (1/sqrt(2*r))*(R2(r-y)(y) - R2(y)(y)))
    S2 = integ(S2 ,2*r-0.5,r*0.5)    
    print(S)
    print(S2)    
    integ1d = 2*(S2+S+tri1)
    print('1d integral', integ1d)

    
    B =  -(3*r+1)/(r*sqrt(r)) + 8*sqrt(1-2*r)/r  + 8*arccos(r/(1-r)) -16*arccos(sqrt(r/(1-r))) / sqrt(r)

    B = 2/sqrt(r)+ (1/pi)*B
    
    
    print('Alt', B)
    return integ1d

def F2(r):
    
    tot_sum = 2 - ((2/pi)*((1/sqrt(r))*arccos((3*r-1)/(1-r)) - arccos(r/(1-r))) + (1/(pi*r))*(-2*sqrt(1-2*r))) - 1/sqrt(r)
    tot_sum = 2 - 3/sqrt(r) + ((2/pi)*((2/sqrt(r))*arcsin(sqrt(r/(1-r))) + arccos(r/(1-r)) + (1/r)*sqrt(1-2*r)))
    tot_sum = 2 - 1/sqrt(r) + (2/pi)*(-(2/sqrt(r))*arccos(sqrt(r/(1-r))) + arccos(r/(1-r)) + (1/r)*sqrt(1-2*r))

    return tot_sum

def test_F3():
    print("F3(r) CDF for r in 1/4 to 1/3")
    r = 0.29
    np,nm=(4000,3999)
    print('F3(%f)='%(r),numerical_CDF_F3(r,np,nm))
    print('Analytic', F3(r))

def test_F2():
    print("F2(r) CDF for r in 1/3 to 1/2")
    r = 0.39
    np,nm=(4000,3999)
    print('F2(%f)='%(r),numerical_CDF_F2(r,np,nm))
    print('Analytic', F2(r))
    
if __name__=='__main__':
#    test_F2_plot()
    test_F2()
#    test_F3()
    
