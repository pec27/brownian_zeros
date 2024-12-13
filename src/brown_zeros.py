"""
Simulating the zeros of brownian motion
Copyright (C) 2024 Peter Creasey
"""

import numpy as np

import pylab as pl

def test_sample_s(t=1.5):
    n = 10000
    u = (0.5+ np.arange(n)) / n
    sin2_ang_u = np.sin(0.5 * np.pi * u)**2
    s = t * 2 * sin2_ang_u / (1+sin2_ang_u)
    
    s_sorted = np.array(sorted(s))
    bin_size = 100
    s_bin = s_sorted[bin_size//2::bin_size]
    s_edges = np.array(list(s_sorted[::bin_size]) + [s_sorted[-1]])
    dens = (bin_size / np.diff(s_edges)) / n

    pl.plot(s_bin, np.log(dens), 'k')
    dens = (1/np.pi)*(t/(2*t-s_bin))*np.sqrt(2/(s_bin * (t-s_bin)))
    pl.plot(s_bin, np.log(dens), 'b')
    pl.show()

def test_dens_r_given_s(s=0.3, t=1.5):

    n = 10000
    dr = t / n
    r = dr * (0.5+np.arange(n))
    dens = np.sqrt((t-s)/(r*t*(2*t-s-r))) * (2*t-s)*0.5 /(2*t-s-r)
    print('total', dens.sum() * dr)

def test_sample_r_given_st(s=0.3, t=1.5):
    n = 10000
    u = (0.5+ np.arange(n)) / n
    u2 = np.square(u)
    r = t * (u2*(2*t-s))/((1+u2)*t - s)
    
    r_sorted = np.array(sorted(r))
    bin_size = 100
    r_bin = r_sorted[bin_size//2::bin_size]
    r_edges = np.array(list(r_sorted[::bin_size]) + [r_sorted[-1]])
    dens = (bin_size / np.diff(r_edges)) / n

    pl.plot(r_bin, np.log(dens), 'k')
    dens = np.sqrt((t-s)/(r_bin*t*(2*t-s-r_bin))) * (2*t-s)*0.5 /(2*t-s-r_bin)    
    pl.plot(r_bin, np.log(dens), 'b')
    pl.show()

def sim(last_sample, threshold, rs):


    zeros = [(0, last_sample)]
    small_intervals = []
    while len(zeros):
        z0,z1 = zeros.pop()
        T = z1 - z0
        if T < threshold:
            print(z0,z1, 'small enough')
            small_intervals.append((z0,z1))

            continue

        u,v = rs.rand(2)
        print(u,v)
        sin2 = np.sin(0.5 * np.pi * u)**2
        s = T * sin2 / (1+sin2)
        v2 = v*v
        
        r = T * (v2*(T-s))/((1+v2)*T - 2*s)
        print('Splitting (%f,%f) to (%f,%f) (%f, %f) (%.2f%% and %.2f%%)'%(z0,z1, z0,z0+s, z1-r,z1, 100*s/T, 100*(T-r)/T))
        zeros.append((z0, z0 + s))
        zeros.append((z1 - r, z1))
    return small_intervals


    
def draw_samples(t_end=100):

    from numpy.random import RandomState
    rs = RandomState(seed=123)
    threshold = 0.001

    lw = 0.1
    c = 'k'
    pl.figure(figsize=(24,5))
    for i in range(5):
        last_sample = t_end #* float(np.square(np.sin(0.5 * np.pi*rs.rand(1))))        
        small_intervals = sim(last_sample, threshold, rs)
        ax = pl.subplot(5,1,i+1)
        for z0,z1 in small_intervals:
            ax.axvline(z0, lw=lw,c=c)
            ax.axvline(z1, lw=lw,c=c)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)        
        pl.xlim(0,t_end)
#        ax.axis('off')
    pl.tight_layout()
#    pl.show()
    pl.savefig('figs/sample.png')
#test_dens_r_given_s(0.1)
#test_sample_r_given_st()
draw_samples()
