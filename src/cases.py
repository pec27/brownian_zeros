"""
Copyright (C) 2024 Peter Creasey
"""

import pylab as pl

pl.rcParams.update({
    "text.usetex": True,
    "font.size": 8,
    "font.family": "serif"#"Helvetica"
})

fontsize=7

fh = 0.013
lc = 'k'
lw = 0.5

def plot_even(n = 3):
    r = (0.4/(2*n) + 0.6/(2*n+1))
    pl.title(r'$F_{2n}(r)$')#,\;r^{-1} \in (2n,2n+1]$')

    pl.plot((0.5+r,0.5),(0.5,0.5+r), lc, lw=lw)
    

    pl.xlim(0.5, 0.5+r*1)
    pl.ylim(0.5, 0.5+r*1)

    v = ([0.5,1-n*r, 0.5+r], [r'$\frac{1}{2}$',r'$1-nr$',r'$\frac{1}{2} + r$'])
    pl.plot((1-n*r,)*2, [0.5,(n+1)*r], lc, lw=lw)
    pl.plot([0.5,(n+1)*r], (1-n*r,)*2, lc, lw=lw)


    # Both n-1
    if True:
        x = 1-n*r + 0.001
        pl.text(x,x+fh,r'$F_{n-1}\left(\frac{r}{\tau_{\mathrm{-}}}\right)$', ha='left',va='bottom',size=fontsize)    
        pl.text(x,x,r'$F_{n-1}\left(\frac{r}{1-\tau_{\mathrm{+}}}\right)$', ha='left',va='bottom',size=fontsize)


    # One n-1, one n
    if True:
        x = 1-n*r + fh*0.1 + 0.01
        y = 0.501#0.3*(1-n*r) + 0.7 * 0.5
        pl.text(x,y+fh,r'$F_{n-1}\left(\frac{r}{\tau_{-}}\right)\times$', ha='left',va='bottom',size=fontsize)        
        pl.text(x,y,r'$F_n\left(\frac{r}{1-\tau_{+}}\right)$', ha='left',va='bottom',size=fontsize)


        x = 0.501
        y = 1-n*r + 0.001
        pl.text(x,y,r'$F_n\left(\frac{r}{\tau_{-}}\right)\times$', ha='left',va='bottom',size=fontsize,rotation=90)    
        pl.text(x+fh,y,r'$F_{n-1}\left(\frac{r}{1-\tau_{+}}\right)$', ha='left',va='bottom',size=fontsize,rotation=90)

    # Both n
    if True:
        x = 0.501#0.5*(1.5-n*r)
        pl.text(x,x+fh,r'$F_n\left(\frac{r}{\tau_{\mathrm{-}}}\right)\times$', ha='left',va='bottom',size=fontsize)            
        pl.text(x,x,r'$F_n\left(\frac{r}{1-\tau_{\mathrm{+}}}\right)$', ha='left',va='bottom',size=fontsize)

    if True:
        x = 0.75 * (0.5+r) + 0.25 * 0.5
        pl.text(x,x,r'$0$', size=fontsize, ha='left', va='bottom')
        
    pl.xticks(v[0],v[1])
    pl.yticks(v[0],v[1])
    pl.xlabel(r'$1-\tau_{-}$')
    pl.ylabel(r'$\tau_{+}$')
    pl.tight_layout()
def plot_odd(n = 3):
    r = (0.91/(2*n+1) + 0.09/(2*n+2))
    pl.title(r'$F_{2n+1}(r)$')#,\;r^{-1} \in (2n+1,2n+2]$')

    pl.plot((0.5+r,0.5),(0.5,0.5+r), lc, lw=lw)
    

    pl.xlim(0.5, 0.5+r*1)
    pl.ylim(0.5, 0.5+r*1)

    v = ([0.5,1-n*r, 0.5+r], [r'$\frac{1}{2}$',r'$1-nr$',r'$\frac{1}{2} + r$'])
    pl.plot((1-n*r,)*2, [0.5,(n+1)*r], lc, lw=lw)
    pl.plot([0.5,(n+1)*r], (1-n*r,)*2, lc, lw=lw)


    # One n-1, one n
    if True:
        x = 1-n*r + 0.001 #0.5*((1-n*r)+0.5*(0.5+r + (n+1)*r))
        y = 0.501#0.3*(1-n*r) + 0.7 * 0.5
        pl.text(x,y+fh,r'$F_{n-1}\left(\frac{r}{\tau_{-}}\right)$', ha='left',va='bottom',size=fontsize)
        pl.text(x,y,r'$F_n\left(\frac{r}{1-\tau_{+}}\right)$', ha='left',va='bottom',size=fontsize)


        x = 0.501
        y = 1-n*r + 0.001
        pl.text(x,y+fh,r'$F_n\left(\frac{r}{\tau_{-}}\right)\times$', ha='left',va='bottom',size=fontsize,rotation=0)    
        pl.text(x,y,r'$F_{n-1}\left(\frac{r}{1-\tau_{+}}\right)$', ha='left',va='bottom',size=fontsize,rotation=0)

    # Both n
    if True:
        x = 0.501#0.5*(1.5-n*r)
        pl.text(x,x+fh,r'$F_n\left(\frac{r}{\tau_{\mathrm{-}}}\right)\times$', ha='left',va='bottom',size=fontsize)    
        pl.text(x,x,r'$F_n\left(\frac{r}{1-\tau_{\mathrm{+}}}\right)$', ha='left',va='bottom',size=fontsize)

    if True:
        x = 0.75 * (0.5+r) + 0.25 * 0.5
        pl.text(x,x,r'$0$', size=fontsize, ha='left', va='bottom')
    pl.xticks(v[0],v[1])
    pl.yticks(v[0],v[1])
    pl.xlabel(r'$1-\tau_{-}$')
    
if __name__=='__main__':
    pl.figure(figsize=(5.8,3),dpi=200)

    # Even
    pl.subplot(1,2,1)
    plot_even()

    pl.subplot(1,2,2)
    plot_odd()
    
    pl.tight_layout()
    pl.savefig('figs/integrand.png', dpi=600)
#    pl.show()
