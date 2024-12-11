import numpy as np
import control as ct
import matplotlib as mpl
import matplotlib.pyplot as plt

s = ct.tf('s')
G = 3*10**4/(s*(s**2+150*s+10**4))

def tidssvar(F):
    L = [Fi*G for Fi in F]
    fig, ax = plt.subplots(2, 1)
    for i in range(len(F)):
        t, er = ct.step_response(30/(s*(1+L[i])), 0.5)  # beräkna rampsvaret enligt krav 1
        ax[0].plot(t, er)
        t, yv = ct.step_response(15*G/(1+L[i]), 0.5)  # beräkna stegsvaret enligt krav 2
        ax[1].plot(t, yv)
        
    fig.suptitle('Slutna systemets tidssvar')
    fig.supxlabel('t')
    plt.subplot(2, 1, 1)
    plt.ylabel('Rampfel $e_r$')
    plt.axhline(0.5, linestyle='--', linewidth=1)
    plt.legend([f'F[{i}]' for i in range(len(F))])
    plt.subplot(2, 1, 2)
    plt.ylabel('Positionsfel $y_v$')
    plt.axhline(0.5, linestyle='--', linewidth=1)
    plt.legend([f'F[{i}]' for i in range(len(F))])
    plt.show()

def Lbode(F):
    L = [Fi*G for Fi in F]
    fig = plt.figure()
    wc = np.zeros(len(F))
    pm = np.zeros(len(F))
    for i in range(len(F)):
        ct.bode_plot(L[i], omega=np.logspace(0, 3))  # plotta Bodediagram för L(s)..
        gmi, pmi, wcgi, wcpi = ct.margin(L[i])      # .. och beräkna fasmarginal och skärfrekvens (krav 3 och 4)
        wc[i] = wcpi
        pm[i] = pmi
    
    fig.suptitle('Bodediagram för kretsöverföringen $L(s)$')
    plt.subplot(2, 1, 1)
    plt.legend([r'$\omega_c=$'+ f'{wc[i]:.0f}' + ' rad/s' for i in range(len(F))])
    plt.subplot(2, 1, 2)
    plt.legend([r'$\varphi_m=$'+ f'{pm[i]:.1f}' + r'$^\circ$' for i in range(len(F))])
    plt.show()

def Tbode(F):
    T = [Fi*G/(1+Fi*G) for Fi in F]
    fig = plt.figure()
    for i in range(len(F)):
        ct.bode_plot(T[i], omega=np.logspace(0, 3))

    fig.suptitle('Bodediagram för slutna systemet $T(s)$')
    plt.show()

Kp1 = 30
F = [Kp1]
tidssvar(F)
Lbode(F)
Tbode(F)

Kp2 = 40/3
a = Kp1/Kp2
tau = 3/40
Flag = a*(1+tau*s)/(1+a*tau*s)
F = [Kp1, Kp2*Flag]
print('Flag = ', Flag)
tidssvar(F)
Lbode(F)
Tbode(F)