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

def dist(d1, d2): 				#position fin colonne
	return (d1+d2)/2

def idist(d1, d2): 				#incertitude sur la position
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
z1a = z(tab1p[63, 0], tab1p[29, 0]) 
k1b = beta(tab1p[63, 1], tab1p[88, 1])
z1b = z(tab1p[63, 0], tab1p[88, 0])
k1c = beta(tab1p[109, 1], tab1p[88, 1])
z1c = z(tab1p[109, 0], tab1p[88, 0])
k1d = beta(tab1p[109, 1], tab1p[127, 1])
z1d = z(tab1p[109, 0], tab1p[127, 0])
k1e = beta(tab1m[15, 1], tab1m[50, 1])
z1e = z(tab1m[15, 0], tab1m[50, 0])
k1f = beta(tab1m[79, 1], tab1m[50, 1])
z1f = z(tab1m[79, 0], tab1m[50, 0])
k1g = beta(tab1m[79, 1], tab1m[104, 1])
z1g = z(tab1m[79, 0], tab1m[104, 0])
k1h = beta(tab1m[124, 1], tab1m[104, 1])
z1h = z(tab1m[124, 0], tab1m[104, 0])
k1j = beta(tab1m[124, 1], tab1m[139, 1])
z1j = z(tab1m[124, 0], tab1m[139, 0])

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

k2a = beta(tab2p[34, 1], tab2p[72, 1])
z2a = z(tab2p[34, 0], tab2p[72, 0])
k2b = beta(tab2p[101, 1], tab2p[72, 1])
z2b = z(tab2p[101, 0], tab2p[72, 0])
k2c = beta(tab2p[101, 1], tab2p[125, 1])
z2c = z(tab2p[101, 0], tab2p[125, 0])
k2d = beta(tab2p[125, 1], tab2p[146, 1])
z2d = z(tab2p[146, 0], tab2p[125, 0])
k2e = beta(tab2p[164, 1], tab2p[146, 1])
z2e = z(tab2p[164, 0], tab2p[146, 0])
k2f = beta(tab2m[11, 1], tab2m[46, 1])
z2f = z(tab2m[11, 0], tab2m[46, 0])
k2g = beta(tab2m[46, 1], tab2m[79, 1])
z2g = z(tab2m[46, 0], tab2m[79,0])
k2h = beta(tab2m[79, 1], tab2m[106, 1])
z2h = z(tab2m[79, 0], tab2m[106,0])
k2i = beta(tab2m[106, 1], tab2m[129, 1])
z2i = z(tab2m[106, 0], tab2m[129,0])
k2j = beta(tab2m[129, 1], tab2m[150, 1])
z2j = z(tab2m[129, 0], tab2m[150,0])  

#Pi = 0.290 & Pr = 0.025
plt.title('3')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab3p[:, 0], tab3p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab3m[:, 0], tab3m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab3g[:, 0], tab3g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

d3 = dist(49.3, 49.2)
id3 = idist(49.3, 49.2)

k3a = beta(tab3p[35, 1], tab3p[75, 1])
z3a = z(tab3p[35, 0], tab3p[75, 0])
k3b = beta(tab3p[75, 1], tab3p[104, 1])
z3b = z(tab3p[75, 0], tab3p[104, 0])
k3c = beta(tab3p[104, 1], tab3p[133, 1])
z3c = z(tab3p[104, 0], tab3p[133, 0])
k3d = beta(tab3p[133, 1], tab3p[158, 1])
z3d = z(tab3p[133, 0], tab3p[158, 0])
k3e = beta(tab3p[158, 1], tab3p[179, 1])
z3e = z(tab3p[158, 0], tab3p[179, 0])
k3f = beta(tab3m[11, 1], tab3m[50, 1])
z3f = z(tab3m[11, 0], tab3m[50, 0])
k3g = beta(tab3m[50, 1], tab3m[86, 1])
z3g = z(tab3m[50, 0], tab3m[86, 0])
k3h = beta(tab3m[86, 1], tab3m[115, 1])
z3h = z(tab3m[86, 0], tab3m[115, 0])
k3i = beta(tab3m[115, 1], tab3m[142, 1])
z3i = z(tab3m[115, 0], tab3m[142, 0])
k3j = beta(tab3m[142, 1], tab3m[163, 1])
z3j = z(tab3m[142, 0], tab3m[163, 0])

#Pi = 0.35 & Pr = 0.021
plt.title('4')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab4p[:, 0], tab4p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab4m[:, 0], tab4m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab4g[:, 0], tab4g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

d4 = dist(45.7, 45.5)
id4 = idist(45.7, 45.5)

k4a = beta(tab4p[30, 1], tab4p[72, 1])
z4a = z(tab4p[30, 0], tab4p[72, 0])
k4b = beta(tab4p[72, 1], tab4p[105, 1])
z4b = z(tab4p[72, 0], tab4p[105, 0])
k4c = beta(tab4p[105, 1], tab4p[133, 1])
z4c = z(tab4p[105, 0], tab4p[133, 0])
k4d = beta(tab4p[133, 1], tab4p[158, 1])
z4d = z(tab4p[133, 0], tab4p[158, 0])
k4e = beta(tab4p[158, 1], tab4p[180, 1])
z4e = z(tab4p[158, 0], tab4p[180, 0])
k4f = beta(tab4m[18, 1], tab4m[53, 1])
z4f = z(tab4m[18, 0], tab4m[53, 0])
k4g = beta(tab4m[53, 1], tab4m[89, 1])
z4g = z(tab4m[53, 0], tab4m[89, 0])
k4h = beta(tab4m[89, 1], tab4m[119, 1])
z4h = z(tab4m[89, 0], tab4m[119, 0])
k4i = beta(tab4m[119, 1], tab4m[147, 1])
z4i = z(tab4m[119, 0], tab4m[147, 0])
k4j = beta(tab4m[147, 1], tab4m[170, 1])
z4j = z(tab4m[147, 0], tab4m[170, 0])

#Pi = 0.40 & Pr = 0.026
plt.title('5')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab5p[:, 0], tab5p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab5m[:, 0], tab5m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab5g[:, 0], tab5g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

d5 = dist(35.3, 34.8)
id5 = idist(35.3, 34.8)

k5a = beta(tab5p[36, 1], tab5p[76, 1])
z5a = z(tab5p[36, 0], tab5p[76, 0])
k5b = beta(tab5p[76, 1], tab5p[114, 1])
z5b = z(tab5p[76, 0], tab5p[114, 0])
k5c = beta(tab5p[114, 1], tab5p[142, 1])
z5c = z(tab5p[114, 0], tab5p[142, 0])
k5d = beta(tab5p[142, 1], tab5p[170, 1])
z5d = z(tab5p[142, 0], tab5p[170, 0])
k5e = beta(tab5p[170, 1], tab5p[195, 1])
z5e = z(tab5p[170, 0], tab5p[195, 0])
k5f = beta(tab5m[15, 1], tab5m[52, 1])
z5f = z(tab5m[15, 0], tab5m[52, 0])
k5g = beta(tab5m[52, 1], tab5m[95, 1])
z5g = z(tab5m[52, 0], tab5m[95, 0])
k5h = beta(tab5m[95, 1], tab5m[129, 1])
z5h = z(tab5m[95, 0], tab5m[129, 0])
k5i = beta(tab5m[129, 1], tab5m[157, 1])
z5i = z(tab5m[129, 0], tab5m[157, 0])
k5j = beta(tab5m[157, 1], tab5m[185, 1])
z5j = z(tab5m[157, 0], tab5m[185, 0])

#Pi = 0.47 & Pr = 0.025
plt.title('6')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab6p[:, 0], tab6p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab6m[:, 0], tab6m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab6g[:, 0], tab6g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

d6 = dist(25.2, 24.4)
id6 = idist(25.2, 24.4)

k6a = beta(tab6p[39, 1], tab6p[84, 1])
z6a = z(tab6p[39, 0], tab6p[84, 0])
k6b = beta(tab6p[84, 1], tab6p[118, 1])
z6b = z(tab6p[84, 0], tab6p[118, 0])
k6c = beta(tab6p[118, 1], tab6p[152, 1])
z6c = z(tab6p[118, 0], tab6p[152, 0])
k6d = beta(tab6p[152, 1], tab6p[185, 1])
z6d = z(tab6p[152, 0], tab6p[185, 0])
k6e = beta(tab6p[185, 1], tab6p[212, 1])
z6e = z(tab6p[185, 0], tab6p[212, 0])
k6f = beta(tab6m[13, 1], tab6m[56, 1])
z6f = z(tab6m[13, 0], tab6m[56, 0])
k6g = beta(tab6m[56, 1], tab6m[98, 1])
z6g = z(tab6m[56, 0], tab6m[98, 0])
k6h = beta(tab6m[98, 1], tab6m[135, 1])
z6h = z(tab6m[98, 0], tab6m[135, 0])
k6i = beta(tab6m[135, 1], tab6m[168, 1])
z6i = z(tab6m[135, 0], tab6m[168, 0])
k6j = beta(tab6m[168, 1], tab6m[198, 1])
z6j = z(tab6m[168, 0], tab6m[198, 0])

#Pi = 0.54 & Pr = 0.028
plt.title('7')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab7m[:, 0], tab7m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab7g[:, 0], tab7g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

d7 = dist(12.8, 13.4)
id7 = idist(12.8, 13.4)

k7a = beta(tab7g[12, 1], tab7g[51, 1])
z7a = z(tab7g[12, 0], tab7g[51, 0])
k7b = beta(tab7g[51, 1], tab7g[97, 1])
z7b = z(tab7g[51, 0], tab7g[97, 0])
k7c = beta(tab7g[97, 1], tab7g[138, 1])
z7c = z(tab7g[97, 0], tab7g[138, 0])
k7d = beta(tab7g[138, 1], tab7g[172, 1])
z7d = z(tab7g[138, 0], tab7g[172, 0])
k7e = beta(tab7g[172, 1], tab7g[206, 1])
z7e = z(tab7g[172, 0], tab7g[206, 0])
k7f = beta(tab7m[15, 1], tab7m[56, 1])
z7f = z(tab7m[15, 0], tab7m[56, 0])
k7g = beta(tab7m[56, 1], tab7m[103, 1])
z7g = z(tab7m[56, 0], tab7m[103, 0])
k7h = beta(tab7m[103, 1], tab7m[140, 1])
z7h = z(tab7m[103, 0], tab7m[140, 0])
k7i = beta(tab7m[140, 1], tab7m[174, 1])
z7i = z(tab7m[140, 0], tab7m[174, 0])
k7j = beta(tab7m[174, 1], tab7m[208, 1])
z7j = z(tab7m[174, 0], tab7m[208, 0])