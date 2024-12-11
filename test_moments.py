import numpy as np

n = 500
u1 = (np.arange(n) + 0.5) / n
a = 0.4
sin2 = np.square(np.sin(u1 * np.pi / 2))

tm = a*sin2 / (1 - a + a*sin2)
m1 = 1-(1-a)**0.5
print('<t-> %.5f, numerically %.5f'%(m1, tm.mean()))

n2 = n*2
u2 = (np.arange(n2) + 0.5) / n2

v = (1-a)*u2*u2
tp = (a + np.multiply.outer(tm, v-1))/(a+np.add.outer(-tm, v))
tp_m1 = a**0.5
print('Average t+ %.5f, numerically %.5f'%(tp_m1, tp.mean()))
print(tm.shape, tp.shape)
av_tp_given_tm_num = tp.mean(axis=1)
av_tp_given_tm = tm - (tm-1)*np.sqrt((a-tm)/(1-a)) * np.arctan(np.sqrt((1-a)/(a-tm)))
print('E[t+|t-]', av_tp_given_tm_num[0], av_tp_given_tm_num[n//2], av_tp_given_tm_num[-1])
print('E[t+|t-]', av_tp_given_tm[0], av_tp_given_tm[n//2], av_tp_given_tm[-1])
tmtp = (np.reshape(tm,(n,1)) * tp).mean(axis=1)
print(tmtp.shape)

m2 = (a**1.5 - (a+2)*np.sqrt(1-a) + 2)/(3)
print('<t-,t+> %.7f numerically %.7f'%(m2, np.inner(np.transpose(tp), tm).mean()/n), tmtp.mean(), (av_tp_given_tm*tm).mean())

# sin integral
ds = a / (n)
s = (0.5 + np.arange(n)) * ds
v = (1/np.pi)*(np.sqrt(s)*np.arcsin(np.sqrt((1-a)/(1-s))))
print('sin integral', v.sum()*ds, (np.sqrt(1-a)+0.5*a*np.sqrt(1-a)+a**1.5-1)/3)
print('scalar integral==<t-^2>', 1 - np.sqrt(1-a) * ((a+2)/2), (tm*tm).mean())

var_tm = 0.5*np.sqrt(1-a)*m1*m1
print('Variance <(t- - <t->)^2>', ((tm-m1)*(tm-m1)).mean(), var_tm)
var_tp = 0.5*np.sqrt(a)*np.square(1-np.sqrt(a))
print('Variance <(t+ - <t+>)^2>', np.square(tp-tp_m1).mean(), var_tp)

covar = m2 - m1*tp_m1 # TODO simplify
print('Covariance <(t+ - <t+>)(t- - <t->)>', covar, ((av_tp_given_tm - tp_m1)*(tm-m1)).mean())

print('Pearson correlation coefficient')
correl = (2/3)*(1-a**0.5 - (1-a)**0.5)/((a*(1-a))**0.25)
print('Correl <(t+ - <t+>)(t- - <t->)>/sqrt(<(t+ - <t+>)^2><(t- - <t->)^2>', correl, ((av_tp_given_tm - tp_m1)*(tm-m1)).mean()/np.sqrt(((tm-m1)*(tm-m1)).mean() *np.square(tp-tp_m1).mean()))