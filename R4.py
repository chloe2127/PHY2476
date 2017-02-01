# -*- coding: utf-8 -*-

###########################################################
# PHY2476 - 4 Densité électronique d'un plasma
# Chloé Lefebvre, Olivier Petit Vincent
# 31 janvier et 7 février 2017
###########################################################

###########################################################
# Importation des modules et des fonctions externes
import numpy as np
import matplotlib.pyplot as plt

###########################################################
# Définition de fonctions locales

def dist(d1, d2):
	return (d1+d2)/2

def idist(d1, d2):
	return np.abs((d1-d2)/2)

def beta(l1, l2):
	l = np.abs(l1-l2)
	return 2*np.pi/l

def z(z1, z2):
	return (z1+z2)/2

###########################################################
# Définition des constantes

# lecture des données txt
diagphase = np.loadtxt('diagramme_de_phase.txt')

tab1p = np.loadtxt('1p.txt')
tab1m = np.loadtxt('1m.txt')
tab1g = np.loadtxt('1g.txt')

tab2p = np.loadtxt('2p.txt')
tab2m = np.loadtxt('2m.txt')
tab2g = np.loadtxt('2g.txt')

tab3p = np.loadtxt('3p.txt')
tab3m = np.loadtxt('3m.txt')
tab3g = np.loadtxt('3g.txt')

tab4p = np.loadtxt('4p.txt')
tab4m = np.loadtxt('4m.txt')
tab4g = np.loadtxt('4g.txt')

tab5p = np.loadtxt('5p.txt')
tab5m = np.loadtxt('5m.txt')
tab5g = np.loadtxt('5g.txt')

tab6p = np.loadtxt('6p.txt')
tab6m = np.loadtxt('6m.txt')
tab6g = np.loadtxt('6g.txt')

tab7m = np.loadtxt('7m.txt')
tab7g = np.loadtxt('7g.txt')

###########################################################
# Programme principal

#diagramme de phase fourni
plt.title("diagramme de phase de l'argon")
plt.xlabel(r'$\beta(m^{-1})$')
plt.ylabel(r'$\omega/\omega_p$')
plt.scatter(diagphase[:, 1], diagphase[:, 0])
plt.show()

#1e5 sur les puissances (en mW)

#Pi = 0.21 & Pr = 0.02
plt.title('1')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab1p[:, 0], tab1p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab1m[:, 0], tab1m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab1g[:, 0], tab1g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

d1 = dist(64.3, 65.2)
id1 = idist(64.3, 65.2)

k1a = beta(tab1p[63, 1],tab1p[29, 1])
z1a = np.mean(tab1p[63, 0], tab1p[29, 0]) 
k1b = beta(tab1p[63, 1], tab1p[88, 1])
z1b = np.mean(tab1p[63, 0], tab1p[88, 0])
k1c = beta(tab1p[109, 1], tap1p[88, 1])
z1c = np.mean(tab1p[109, 0], tab1p[88, 0])
k1d = beta(tab1p[109, 1], tap1p[127, 1])
z1d = np.mean(tab1p[109, 0], tab1p[127, 0])
k1e = beta(tab1m[15, 1], tap1m[50, 1])
z1e = np.mean(tab1m[15, 0], tab1m[50, 0])
k1f = beta(tab1m[79, 1], tap1m[50, 1])
z1f = np.mean(tab1m[79, 0], tab1m[50, 0])
k1g = beta(tab1m[79, 1], tap1m[104, 1])
z1g = np.mean(tab1m[79, 0], tab1m[104, 0])
k1h = beta(tab1m[124, 1], tap1m[104, 1])
z1h = np.mean(tab1m[124, 0], tab1m[104, 0])
k1j = beta(tab1m[124, 1], tap1m[139, 1])
z1j = np.mean(tab1m[124, 0], tab1m[139, 0])

#Pi = 0.245 & Pr = 0.02
plt.title('2')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab2p[:, 0], tab2p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab2m[:, 0], tab2m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab2g[:, 0], tab2g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

d2 = dist(55.3, 55.5)
id2 = idist(55.3, 55.5)

# k2a = beta(tab2p[], tab2p[])
# z2a = np.mean()
# k2b = 
# z2b = 
# k2c = 
# z2c =
# k2d = 
# z2d = 
# k2e = 
# z2e = 
# k2f = 
# z2f = 
# k2g = 
# z2g = 
# k2h = 
# z2h = 
# k2i = 
# z2i = 
# k2j = 
# z2j =  

#Pi = 0.290 & Pr = 0.025
plt.title('3')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab3p[:, 0], tab3p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab3m[:, 0], tab3m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab3g[:, 0], tab3g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

#Pi = 0.35 & Pr = 0.021
plt.title('4')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab4p[:, 0], tab4p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab4m[:, 0], tab4m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab4g[:, 0], tab4g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

#Pi = 0.40 & Pr = 0.026
plt.title('5')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab5p[:, 0], tab5p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab5m[:, 0], tab5m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab5g[:, 0], tab5g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

#Pi = 0.47 & Pr = 0.025
plt.title('6')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab6p[:, 0], tab6p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab6m[:, 0], tab6m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab6g[:, 0], tab6g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

#Pi = 0.54 & Pr = 0.028
plt.title('7')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab7m[:, 0], tab7m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab7g[:, 0], tab7g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()