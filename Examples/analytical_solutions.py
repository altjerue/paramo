import numpy as np
from scipy.special import iv
import matplotlib.pyplot as plt
import Eduviges.constants as acons

# solution from Fokker-Planck Equations of Stochastic Acceleration: Green's Functions and Boundary Conditions by Park and Petrosian 1995

def eq_59_Park1995(t,g):
    D_cof = 1e0

    x = np.sqrt(g**2 -1)#*acons.cLight*acons.me#np.logspace(-3,3,100) #g-1
    q = 3
    a = 1 # a = fancya/fancyd
    x0 = np.sqrt((100**2)- 1)
     #at t=30 the function is scaled p by a factor of 1e11
    alpha = (2-q)/2
    theta=1
    t_esc = 1#fancy F #theta = 1/(D t_esc)
    tau = D_cof * t
    gfa =(1/np.abs(alpha))*(1/(2*tau))
    gfb =(x**((1-q +a)/2))*(x0**((1-q -a)/2))
    etta = (q-1 + a )/(2*alpha)
    order = np.abs(etta)
    gfc = iv(order, (x**alpha)*(x0**alpha)/(2*tau*(alpha**2)))
    gfd = np.exp(-((x**(2*alpha)) + (x0**(2*alpha)))/(4*(alpha**2)*tau))
    gfe = np.exp(-theta*tau)
    Gf = gfa*gfb*gfc*gfd*gfe
    try:
        iter(Gf)
        for i in range(len(x)):
            if(np.isinf(Gf[i]) or np.isnan(Gf[i])):
                Gf[i] = 0e0
    except:
        if(np.isinf(Gf) or np.isnan(Gf)):
            Gf = 0e0
    # Gf=n0

    return Gf*4*np.pi

# x = np.logspace(0,4,100)
# y = eq_59_Park1995(0.003,x)
# fig, ax = plt.subplots()
# ax.set_yscale('log')
# ax.set_xscale('log')
# ax.set_ylim([1e-8,1e2])
# ax.set_xlim([1e-4,1e4])
# pl = ax.plot(x, y)
# ax.legend()
# plt.show()
