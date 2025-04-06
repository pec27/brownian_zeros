"""
Copyright (C) 2024 Peter Creasey
"""
import numpy as np

n = 500
u1 = (np.arange(n) + 0.5) / n
a = 0.4
sin2 = np.square(np.sin(u1 * np.pi / 2))

tp = a / (a + (1-a)*sin2)
tp_m1 = a**0.5
print('<t+> %.5f, numerically %.5f'%(tp_m1, tp.mean()))

n2 = n*2
u2 = (np.arange(n2) + 0.5) / n2

v = a*u2*u2
tm = (np.multiply.outer(tp, v))/(np.add.outer(tp, v)-a)
tm_m1 = 1-(1-a)**0.5
print('<t-> %.5f, numerically %.5f'%(tm_m1, tm.mean()))
print(tm.shape, tp.shape)
av_tm_given_tp_num = tm.mean(axis=1)
av_tm_given_tp = tp - tp*np.sqrt((tp-a)/a) * np.arctan(np.sqrt((a)/(tp-a)))
print('E[t-|t+]', av_tm_given_tp_num[0], av_tm_given_tp_num[n//2], av_tm_given_tp_num[-1])
print('E[t-|t+]', av_tm_given_tp[0], av_tm_given_tp[n//2], av_tm_given_tp[-1])
tmtp = (np.reshape(tp,(n,1)) * tm).mean(axis=1)
print(tmtp.shape)

m2 = (a**1.5 - (a+2)*np.sqrt(1-a) + 2)/(3)
print('<t-,t+> %.7f numerically %.7f'%(m2, np.inner(np.transpose(tm), tp).mean()/n), tmtp.mean(), (av_tm_given_tp*tp).mean())
print('Median t-', a/(2-a), np.median(tm))
print('Median t+', 2*a/(1+a), np.median(tp.ravel()))

# sin integral
ds = a / n
s = (0.5 + np.arange(n)) * ds
v = (1/np.pi)*(np.sqrt(s)*np.arcsin(np.sqrt((1-a)/(1-s))))
print('sin integral', v.sum()*ds, (np.sqrt(1-a)+0.5*a*np.sqrt(1-a)+a**1.5-1)/3)
print('scalar integral==<t-^2>', 1 - np.sqrt(1-a) * ((a+2)/2), (tm*tm).mean())

var_tm = 0.5*np.sqrt(1-a)*np.square(tm_m1)
print('Variance <(t- - <t->)^2>', ((tm-tm_m1)*(tm-tm_m1)).mean(), var_tm)
var_tp = 0.5*np.sqrt(a)*np.square(1-np.sqrt(a))
print('Variance <(t+ - <t+>)^2>', np.square(tp-tp_m1).mean(), var_tp)

covar = m2 - tm_m1*tp_m1 # TODO simplify
print('Covariance <(t+ - <t+>)(t- - <t->)>', covar, ((av_tm_given_tp - tm_m1)*(tp-tp_m1)).mean())

print('Pearson correlation coefficient')
correl = (2/3)*(1-a**0.5 - (1-a)**0.5)/((a*(1-a))**0.25)
print('Correl <(t+ - <t+>)(t- - <t->)>/sqrt(<(t+ - <t+>)^2><(t- - <t->)^2>', correl, ((av_tm_given_tp - tm_m1)*(tp-tp_m1)).mean()/np.sqrt(((tm-tm_m1)*(tm-tm_m1)).mean() *np.square(tp-tp_m1).mean()))


print('Correl a=0.5, corr= (8**0.5 - 4)/3', (8**0.5-4)/3)
