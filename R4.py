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

def dist(d1, d2): 				#position fin colonne (en m)
	return (d1+d2)/200

def idist(d1, d2): 				#incertitude sur la position (en m)
	return np.abs((d1-d2)/200)

def beta(l1, l2): 				#retourne k (beta) le vecteur d'onde (en m-1)
	l = np.abs(l1-l2)*(200)
	return 2*np.pi/l

def z(z1, z2): 					#position antenne (en m)
	return (z1+z2)/200

def densite(x):					#Le gros polynome pour 1/sqrt(n_e)	x = k
	rapport = -1E-12*x**6 + 7E-10*x**5 - 2E-07*x**4 + 2E-05*x**3 - 0.0012*x**2 + 0.0421*x - 0.3292	
	densite = (omega**2)*me*epsilon0/((rapport**2)*(e**2))
	return densite

###########################################################
# Définition des constantes

#constantes pour fct densité
omega = 2*np.pi*600e6			#2pi*f (600 MHz)
me = 9.10938356e-31 			#masse au repos d'un électron (kg)
epsilon0 = 8.8518782e-12		#permittivité du vide (m-3kg-1s4A2)
e = 1.60217662e-19				#charge électronique (C)

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

tab163_p = np.loadtxt('163_p.txt')
tab163_m = np.loadtxt('163_m.txt')
tab163_g = np.loadtxt('163_g.txt')

tab261_p = np.loadtxt('261_p.txt')
tab261_m = np.loadtxt('261_m.txt')
tab261_g = np.loadtxt('261_g.txt')

tab363_p = np.loadtxt('363_p.txt')
tab363_m = np.loadtxt('363_m.txt')
tab363_g = np.loadtxt('363_g.txt')

#tab455_p = np.loadtxt('455_p.txt')
tab455_m = np.loadtxt('455_m.txt')
tab455_g = np.loadtxt('455_g.txt')

tab573_p = np.loadtxt('573_p.txt')
tab573_m = np.loadtxt('573_m.txt')
tab573_g = np.loadtxt('573_g.txt')

tab665_p = np.loadtxt('665_p.txt')
tab665_m = np.loadtxt('665_m.txt')
tab665_g = np.loadtxt('665_g.txt')

tab766_p = np.loadtxt('766_p.txt')
tab766_m = np.loadtxt('766_m.txt')
tab766_g = np.loadtxt('766_g.txt')

tab880_p = np.loadtxt('880_p.txt')
tab880_m = np.loadtxt('880_m.txt')
tab880_g = np.loadtxt('880_g.txt')

tab1024_p = np.loadtxt('1024_p.txt')
tab1024_m = np.loadtxt('1024_m.txt')
tab1024_g = np.loadtxt('1024_g.txt')

###########################################################
# Programme principal

#diagramme de phase fourni
plt.title("diagramme de phase de l'argon")
plt.xlabel(r'$\beta(m^{-1})$')
plt.ylabel(r'$\omega/\omega_p$')
plt.scatter(diagphase[:, 1], diagphase[:, 0])
plt.show()

#1e5 sur les puissances (en mW)

#A PRESSION CONSTANTE (50 mTorrs)

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

k1a = beta(tab1p[29, 0],tab1p[62, 0])
z1a = z(tab1p[29, 0], tab1p[62, 0]) 
k1b = beta(tab1p[63, 0], tab1p[88, 0])
z1b = z(tab1p[63, 0], tab1p[88, 0])
k1c = beta(tab1p[109, 0], tab1p[88, 0])
z1c = z(tab1p[109, 0], tab1p[88, 0])
k1d = beta(tab1p[109, 0], tab1p[127, 0])
z1d = z(tab1p[109, 0], tab1p[127, 0])
k1e = beta(tab1m[15, 0], tab1m[50, 0])
z1e = z(tab1m[15, 0], tab1m[50, 0])
k1f = beta(tab1m[79, 0], tab1m[50, 0])
z1f = z(tab1m[79, 0], tab1m[50, 0])
k1g = beta(tab1m[79, 0], tab1m[104, 0])
z1g = z(tab1m[79, 0], tab1m[104, 0])
k1h = beta(tab1m[124, 0], tab1m[104, 0])
z1h = z(tab1m[124, 0], tab1m[104, 0])
k1j = beta(tab1m[124, 0], tab1m[139, 0])
z1j = z(tab1m[124, 0], tab1m[139, 0])

tabk1 = np.zeros((9,2))
tabk1[0,:] = k1a, z1a
tabk1[1,:] = k1b, z1b
tabk1[2,:] = k1c, z1c
tabk1[3,:] = k1d, z1d
tabk1[4,:] = k1e, z1e
tabk1[5,:] = k1f, z1f
tabk1[6,:] = k1g, z1g
tabk1[7,:] = k1h, z1h
tabk1[8,:] = k1j, z1j
print(tabk1)

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

k2a = beta(tab2p[34, 0], tab2p[72, 0])
z2a = z(tab2p[34, 0], tab2p[72, 0])
k2b = beta(tab2p[101, 0], tab2p[72, 0])
z2b = z(tab2p[101, 0], tab2p[72, 0])
k2c = beta(tab2p[101, 0], tab2p[125, 0])
z2c = z(tab2p[101, 0], tab2p[125, 0])
k2d = beta(tab2p[125, 0], tab2p[146, 0])
z2d = z(tab2p[146, 0], tab2p[125, 0])
k2e = beta(tab2p[164, 0], tab2p[146, 0])
z2e = z(tab2p[164, 0], tab2p[146, 0])
k2f = beta(tab2m[11, 0], tab2m[46, 0])
z2f = z(tab2m[11, 0], tab2m[46, 0])
k2g = beta(tab2m[46, 0], tab2m[79, 0])
z2g = z(tab2m[46, 0], tab2m[79,0])
k2h = beta(tab2m[79, 0], tab2m[106, 0])
z2h = z(tab2m[79, 0], tab2m[106,0])
k2i = beta(tab2m[106, 0], tab2m[129, 0])
z2i = z(tab2m[106, 0], tab2m[129,0])
k2j = beta(tab2m[129, 0], tab2m[150, 0])
z2j = z(tab2m[129, 0], tab2m[150,0])  

tabk2 = np.zeros((9,2))
tabk2[0,:] = k2a, z2a
tabk2[1,:] = k2b, z2b
tabk2[2,:] = k2c, z2c
tabk2[3,:] = k2d, z2d
tabk2[4,:] = k2e, z2e
tabk2[5,:] = k2f, z2f
tabk2[6,:] = k2g, z2g
tabk2[7,:] = k2h, z2h
tabk2[8,:] = k2j, z2j

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

k3a = beta(tab3p[35, 0], tab3p[75, 0])
z3a = z(tab3p[35, 0], tab3p[75, 0])
k3b = beta(tab3p[75, 0], tab3p[104, 0])
z3b = z(tab3p[75, 0], tab3p[104, 0])
k3c = beta(tab3p[104, 0], tab3p[133, 0])
z3c = z(tab3p[104, 0], tab3p[133, 0])
k3d = beta(tab3p[133, 0], tab3p[158, 0])
z3d = z(tab3p[133, 0], tab3p[158, 0])
k3e = beta(tab3p[158, 0], tab3p[179, 0])
z3e = z(tab3p[158, 0], tab3p[179, 0])
k3f = beta(tab3m[11, 0], tab3m[50, 0])
z3f = z(tab3m[11, 0], tab3m[50, 0])
k3g = beta(tab3m[50, 0], tab3m[86, 0])
z3g = z(tab3m[50, 0], tab3m[86, 0])
k3h = beta(tab3m[86, 0], tab3m[115, 0])
z3h = z(tab3m[86, 0], tab3m[115, 0])
k3i = beta(tab3m[115, 0], tab3m[142, 0])
z3i = z(tab3m[115, 0], tab3m[142, 0])
k3j = beta(tab3m[142, 0], tab3m[163, 0])
z3j = z(tab3m[142, 0], tab3m[163, 0])

tabk3 = np.zeros((9,2))
tabk3[0,:] = k3a, z3a
tabk3[1,:] = k3b, z3b
tabk3[2,:] = k3c, z3c
tabk3[3,:] = k3d, z3d
tabk3[4,:] = k3e, z3e
tabk3[5,:] = k3f, z3f
tabk3[6,:] = k3g, z3g
tabk3[7,:] = k3h, z3h
tabk3[8,:] = k3j, z3j

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

k4a = beta(tab4p[30, 0], tab4p[72, 0])
z4a = z(tab4p[30, 0], tab4p[72, 0])
k4b = beta(tab4p[72, 0], tab4p[105, 0])
z4b = z(tab4p[72, 0], tab4p[105, 0])
k4c = beta(tab4p[105, 0], tab4p[133, 0])
z4c = z(tab4p[105, 0], tab4p[133, 0])
k4d = beta(tab4p[133, 0], tab4p[158, 0])
z4d = z(tab4p[133, 0], tab4p[158, 0])
k4e = beta(tab4p[158, 0], tab4p[180, 0])
z4e = z(tab4p[158, 0], tab4p[180, 0])
k4f = beta(tab4m[18, 0], tab4m[53, 0])
z4f = z(tab4m[18, 0], tab4m[53, 0])
k4g = beta(tab4m[53, 0], tab4m[89, 0])
z4g = z(tab4m[53, 0], tab4m[89, 0])
k4h = beta(tab4m[89, 0], tab4m[119, 0])
z4h = z(tab4m[89, 0], tab4m[119, 0])
k4i = beta(tab4m[119, 0], tab4m[147, 0])
z4i = z(tab4m[119, 0], tab4m[147, 0])
k4j = beta(tab4m[147, 0], tab4m[170, 0])
z4j = z(tab4m[147, 0], tab4m[170, 0])

tabk4 = np.zeros((9,2))
tabk4[0,:] = k4a, z4a
tabk4[1,:] = k4b, z4b
tabk4[2,:] = k4c, z4c
tabk4[3,:] = k4d, z4d
tabk4[4,:] = k4e, z4e
tabk4[5,:] = k4f, z4f
tabk4[6,:] = k4g, z4g
tabk4[7,:] = k4h, z4h
tabk4[8,:] = k4j, z4j

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

k5a = beta(tab5p[36, 0], tab5p[76, 0])
z5a = z(tab5p[36, 0], tab5p[76, 0])
k5b = beta(tab5p[76, 0], tab5p[114, 0])
z5b = z(tab5p[76, 0], tab5p[114, 0])
k5c = beta(tab5p[114, 0], tab5p[142, 0])
z5c = z(tab5p[114, 0], tab5p[142, 0])
k5d = beta(tab5p[142, 0], tab5p[170, 0])
z5d = z(tab5p[142, 0], tab5p[170, 0])
k5e = beta(tab5p[170, 0], tab5p[195, 0])
z5e = z(tab5p[170, 0], tab5p[195, 0])
k5f = beta(tab5m[15, 0], tab5m[52, 0])
z5f = z(tab5m[15, 0], tab5m[52, 0])
k5g = beta(tab5m[52, 0], tab5m[95, 0])
z5g = z(tab5m[52, 0], tab5m[95, 0])
k5h = beta(tab5m[95, 0], tab5m[129, 0])
z5h = z(tab5m[95, 0], tab5m[129, 0])
k5i = beta(tab5m[129, 0], tab5m[157, 0])
z5i = z(tab5m[129, 0], tab5m[157, 0])
k5j = beta(tab5m[157, 0], tab5m[185, 0])
z5j = z(tab5m[157, 0], tab5m[185, 0])

tabk5 = np.zeros((9,2))
tabk5[0,:] = k5a, z5a
tabk5[1,:] = k5b, z5b
tabk5[2,:] = k5c, z5c
tabk5[3,:] = k5d, z5d
tabk5[4,:] = k5e, z5e
tabk5[5,:] = k5f, z5f
tabk5[6,:] = k5g, z5g
tabk5[7,:] = k5h, z5h
tabk5[8,:] = k5j, z5j

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

k6a = beta(tab6p[39, 0], tab6p[84, 0])
z6a = z(tab6p[39, 0], tab6p[84, 0])
k6b = beta(tab6p[84, 0], tab6p[118, 0])
z6b = z(tab6p[84, 0], tab6p[118, 0])
k6c = beta(tab6p[118, 0], tab6p[152, 0])
z6c = z(tab6p[118, 0], tab6p[152, 0])
k6d = beta(tab6p[152, 0], tab6p[185, 0])
z6d = z(tab6p[152, 0], tab6p[185, 0])
k6e = beta(tab6p[185, 0], tab6p[212, 0])
z6e = z(tab6p[185, 0], tab6p[212, 0])
k6f = beta(tab6m[13, 0], tab6m[56, 0])
z6f = z(tab6m[13, 0], tab6m[56, 0])
k6g = beta(tab6m[56, 0], tab6m[98, 0])
z6g = z(tab6m[56, 0], tab6m[98, 0])
k6h = beta(tab6m[98, 0], tab6m[135, 0])
z6h = z(tab6m[98, 0], tab6m[135, 0])
k6i = beta(tab6m[135, 0], tab6m[168, 0])
z6i = z(tab6m[135, 0], tab6m[168, 0])
k6j = beta(tab6m[168, 0], tab6m[198, 0])
z6j = z(tab6m[168, 0], tab6m[198, 0])

tabk6 = np.zeros((9,2))
tabk6[0,:] = k6a, z6a
tabk6[1,:] = k6b, z6b
tabk6[2,:] = k6c, z6c
tabk6[3,:] = k6d, z6d
tabk6[4,:] = k6e, z6e
tabk6[5,:] = k6f, z6f
tabk6[6,:] = k6g, z6g
tabk6[7,:] = k6h, z6h
tabk6[8,:] = k6j, z6j

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

k7a = beta(tab7g[12, 0], tab7g[51, 0])
z7a = z(tab7g[12, 0], tab7g[51, 0])
k7b = beta(tab7g[51, 0], tab7g[97, 0])
z7b = z(tab7g[51, 0], tab7g[97, 0])
k7c = beta(tab7g[97, 0], tab7g[138, 0])
z7c = z(tab7g[97, 0], tab7g[138, 0])
k7d = beta(tab7g[138, 0], tab7g[172, 0])
z7d = z(tab7g[138, 0], tab7g[172, 0])
k7e = beta(tab7g[172, 0], tab7g[206, 0])
z7e = z(tab7g[172, 0], tab7g[206, 0])
k7f = beta(tab7m[15, 0], tab7m[56, 0])
z7f = z(tab7m[15, 0], tab7m[56, 0])
k7g = beta(tab7m[56, 0], tab7m[103, 0])
z7g = z(tab7m[56, 0], tab7m[103, 0])
k7h = beta(tab7m[103, 0], tab7m[140, 0])
z7h = z(tab7m[103, 0], tab7m[140, 0])
k7i = beta(tab7m[140, 0], tab7m[174, 0])
z7i = z(tab7m[140, 0], tab7m[174, 0])
k7j = beta(tab7m[174, 0], tab7m[208, 0])
z7j = z(tab7m[174, 0], tab7m[208, 0])

tabk7 = np.zeros((9,2))
tabk7[0,:] = k7a, z7a
tabk7[1,:] = k7b, z7b
tabk7[2,:] = k7c, z7c
tabk7[3,:] = k7d, z7d
tabk7[4,:] = k7e, z7e
tabk7[5,:] = k7f, z7f
tabk7[6,:] = k7g, z7g
tabk7[7,:] = k7h, z7h
tabk7[8,:] = k7j, z7j

#calcul des densités électroniques

plt.title('densité électronique')
plt.xlabel('position (cm)')
plt.ylabel('densité électronique (???)')
plt.scatter(tabk1[:, 1], densite(tabk1[:, 0]), color = 'c', marker = '.', label='1')
#plt.scatter(tabk2[:,1], densite(tabk2[:, 0]), color = 'b', marker = '.', label='2')
#plt.scatter(tabk3[:,1], densite(tabk3[:, 0]), color = 'r', marker = '.', label='3')
plt.legend(loc='best')
plt.show()

######################################################################################
#1e5 sur les puissances (en mW)

#A PUISSANCE 'CONSTANTE' 

#Pi = 0.56 & Pr = 0.04
plt.title('163')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab163_p[:, 0], tab163_p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab163_m[:, 0], tab163_m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab163_g[:, 0], tab163_g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

d163 = dist(24, 23.9)
id163 = idist(24, 23.9)

tabk163 = np.zeros((9,2))
tabk163[0, 0] = beta(tab163_p[35, 0],tab163_p[91, 0])
tabk163[0, 1] = z(tab163_p[35, 0], tab163_p[91, 0]) 
tabk163[1, 0] = beta(tab163_p[91, 0], tab163_p[143, 0])
tabk163[1, 1] = z(tab163_p[91, 0], tab163_p[143, 0])
tabk163[2, 0] = beta(tab163_p[143, 0], tab163_p[187, 0])
tabk163[2, 1] = z(tab163_p[143, 0], tab163_p[187, 0])
tabk163[3, 0] = beta(tab163_p[187, 0], tab163_p[220, 0])
tabk163[3, 1] = z(tab163_p[187, 0], tab163_p[220, 0])
tabk163[4, 0] = beta(tab163_m[66, 0], tab163_m[119, 0])
tabk163[4, 1] = z(tab163_m[66, 0], tab163_m[119, 0])
tabk163[5, 0] = beta(tab163_m[119, 0], tab163_m[167, 0])
tabk163[5, 1] = z(tab163_m[119, 0], tab163_m[167, 0])
tabk163[6, 0] = beta(tab163_m[167, 0], tab163_m[207, 0])
tabk163[6, 1] = z(tab163_m[167, 0], tab163_m[207, 0])
tabk163[7, 0] = beta(tab163_m[207, 0], tab163_m[245, 0])
tabk163[7, 1] = z(tab163_m[207, 0], tab163_m[245, 0])
tabk163[8, 0] = beta(tab163_m[245, 0], tab163_m[271, 0])
tabk163[8, 1] = z(tab163_m[245, 0], tab163_m[271, 0])

#Pi = 0.73 & Pr = 0.02
plt.title('261')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab261_p[:, 0], tab261_p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab261_m[:, 0], tab261_m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab261_g[:, 0], tab261_g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

d261 = dist(17.7, 17.2)
id261 = idist(17.7, 17.2)

tabk261 = np.zeros((9,2))
tabk261[0, 0] = beta(tab261_p[36, 0],tab261_p[106, 0])
tabk261[0, 1] = z(tab261_p[36, 0], tab261_p[106, 0]) 
tabk261[1, 0] = beta(tab261_p[106, 0], tab261_p[169, 0])
tabk261[1, 1] = z(tab261_p[106, 0], tab261_p[169, 0])
tabk261[2, 0] = beta(tab261_p[169, 0], tab261_p[216, 0])
tabk261[2, 1] = z(tab261_p[169, 0], tab261_p[216, 0])
tabk261[3, 0] = beta(tab261_p[216, 0], tab261_p[259, 0])
tabk261[3, 1] = z(tab261_p[216, 0], tab261_p[259, 0])
tabk261[4, 0] = beta(tab261_m[16, 0], tab261_m[79, 0])
tabk261[4, 1] = z(tab261_m[16, 0], tab261_m[79, 0])
tabk261[5, 0] = beta(tab261_m[79, 0], tab261_m[141, 0])
tabk261[5, 1] = z(tab261_m[79, 0], tab261_m[141, 0])
tabk261[6, 0] = beta(tab261_m[141, 0], tab261_m[197, 0])
tabk261[6, 1] = z(tab261_m[141, 0], tab261_m[197, 0])
tabk261[7, 0] = beta(tab261_m[197, 0], tab261_m[241, 0])
tabk261[7, 1] = z(tab261_m[197, 0], tab261_m[241, 0])
tabk261[8, 0] = beta(tab261_m[241, 0], tab261_m[282, 0])
tabk261[8, 1] = z(tab261_m[241, 0], tab261_m[282, 0])

#Pi = 0.72 & Pr = 0.01
plt.title('363')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab363_p[:, 0], tab363_p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab363_m[:, 0], tab363_m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab363_g[:, 0], tab363_g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

d363 = dist(49.3, 49.2)
id363 = idist(49.3, 49.2)

tabk363 = np.zeros((9,2))
tabk363[0, 0] = beta(tab363_p[42, 0],tab363_p[107, 0])
tabk363[0, 1] = z(tab363_p[42, 0], tab363_p[107, 0]) 
tabk363[1, 0] = beta(tab363_p[107, 0], tab363_p[172, 0])
tabk363[1, 1] = z(tab363_p[107, 0], tab363_p[172, 0])
tabk363[2, 0] = beta(tab363_p[172, 0], tab363_p[221, 0])
tabk363[2, 1] = z(tab363_p[172, 0], tab363_p[221, 0])
tabk363[3, 0] = beta(tab363_p[221, 0], tab363_p[262, 0])
tabk363[3, 1] = z(tab363_p[221, 0], tab363_p[262, 0])
tabk363[4, 0] = beta(tab363_m[13, 0], tab363_m[70, 0])
tabk363[4, 1] = z(tab363_m[13, 0], tab363_m[70, 0])
tabk363[5, 0] = beta(tab363_m[70, 0], tab363_m[132, 0])
tabk363[5, 1] = z(tab363_m[70, 0], tab363_m[132, 0])
tabk363[6, 0] = beta(tab363_m[132, 0], tab363_m[191, 0])
tabk363[6, 1] = z(tab363_m[132, 0], tab363_m[191, 0])
tabk363[7, 0] = beta(tab363_m[191, 0], tab363_m[236, 0])
tabk363[7, 1] = z(tab363_m[191, 0], tab363_m[236, 0])
tabk363[8, 0] = beta(tab363_m[236, 0], tab363_m[272, 0])
tabk363[8, 1] = z(tab363_m[236, 0], tab363_m[272, 0])

#Pi = 0.69 & Pr = 0
plt.title('455')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
#plt.scatter(tab455_p[:, 0], tab455_p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab455_m[:, 0], tab455_m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab455_g[:, 0], tab455_g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

d455 = dist(21.3, 20.7)
id455 = idist(21.3, 20.7)

tabk455 = np.zeros((9,2))
tabk455[0, 0] = beta(tab455_g[63, 0],tab455_g[128, 0])
tabk455[0, 1] = z(tab455_g[63, 0], tab455_g[128, 0]) 
tabk455[1, 0] = beta(tab455_g[128, 0], tab455_g[196, 0])
tabk455[1, 1] = z(tab455_g[128, 0], tab455_g[196, 0])
tabk455[2, 0] = beta(tab455_g[196, 0], tab455_g[245, 0])
tabk455[2, 1] = z(tab455_g[196, 0], tab455_g[245, 0])
tabk455[3, 0] = beta(tab455_g[245, 0], tab455_g[287, 0])
tabk455[3, 1] = z(tab455_g[245, 0], tab455_g[287, 0])
tabk455[4, 0] = beta(tab455_m[16, 0], tab455_m[73, 0])
tabk455[4, 1] = z(tab455_m[16, 0], tab455_m[73, 0])
tabk455[5, 0] = beta(tab455_m[73, 0], tab455_m[141, 0])
tabk455[5, 1] = z(tab455_m[73, 0], tab455_m[141, 0])
tabk455[6, 0] = beta(tab455_m[141, 0], tab455_m[206, 0])
tabk455[6, 1] = z(tab455_m[141, 0], tab455_m[206, 0])
tabk455[7, 0] = beta(tab455_m[206, 0], tab455_m[256, 0])
tabk455[7, 1] = z(tab455_m[206, 0], tab455_m[256, 0])
tabk455[8, 0] = beta(tab455_m[256, 0], tab455_m[295, 0])
tabk455[8, 1] = z(tab455_m[256, 0], tab455_m[295, 0])

#Pi = 0.70 & Pr = 0.01
plt.title('573')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab573_p[:, 0], tab573_p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab573_m[:, 0], tab573_m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab573_g[:, 0], tab573_g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

d573 = dist(17.7, 16.8)
id573 = idist(17.7, 16.8)

tabk573 = np.zeros((9,2))
tabk573[0, 0] = beta(tab573_p[46, 0],tab573_p[116, 0])
tabk573[0, 1] = z(tab573_p[46, 0], tab573_p[116, 0]) 
tabk573[1, 0] = beta(tab573_p[116, 0], tab573_p[190, 0])
tabk573[1, 1] = z(tab573_p[116, 0], tab573_p[190, 0])
tabk573[2, 0] = beta(tab573_p[190, 0], tab573_p[249, 0])
tabk573[2, 1] = z(tab573_p[190, 0], tab573_p[249, 0])
tabk573[3, 0] = beta(tab573_p[249, 0], tab573_p[298, 0])
tabk573[3, 1] = z(tab573_p[249, 0], tab573_p[298, 0])
tabk573[4, 0] = beta(tab573_m[19, 0], tab573_m[75, 0])
tabk573[4, 1] = z(tab573_m[19, 0], tab573_m[75, 0])
tabk573[5, 0] = beta(tab573_m[75, 0], tab573_m[146, 0])
tabk573[5, 1] = z(tab573_m[75, 0], tab573_m[146, 0])
tabk573[6, 0] = beta(tab573_m[146, 0], tab573_m[214, 0])
tabk573[6, 1] = z(tab573_m[146, 0], tab573_m[214, 0])
tabk573[7, 0] = beta(tab573_m[214, 0], tab573_m[270, 0])
tabk573[7, 1] = z(tab573_m[214, 0], tab573_m[270, 0])
tabk573[8, 0] = beta(tab573_m[270, 0], tab573_m[306, 0])
tabk573[8, 1] = z(tab573_m[270, 0], tab573_m[306, 0])

#Pi = 0.72 & Pr = 0.04 
plt.title('665')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab665_p[:, 0], tab665_p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab665_m[:, 0], tab665_m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab665_g[:, 0], tab665_g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

d665 = dist(15.3, 15.5)
id665 = idist(15.3, 15.5)

tabk665 = np.zeros((9,2))
tabk665[0, 0] = beta(tab665_p[41, 0],tab665_p[113, 0])
tabk665[0, 1] = z(tab665_p[41, 0], tab665_p[113, 0]) 
tabk665[1, 0] = beta(tab665_p[113, 0], tab665_p[186, 0])
tabk665[1, 1] = z(tab665_p[113, 0], tab665_p[186, 0])
tabk665[2, 0] = beta(tab665_p[186, 0], tab665_p[243, 0])
tabk665[2, 1] = z(tab665_p[186, 0], tab665_p[243, 0])
tabk665[3, 0] = beta(tab665_p[243, 0], tab665_p[296, 0])
tabk665[3, 1] = z(tab665_p[243, 0], tab665_p[296, 0])
tabk665[4, 0] = beta(tab665_m[15, 0], tab665_m[75, 0])
tabk665[4, 1] = z(tab665_m[15, 0], tab665_m[75, 0])
tabk665[5, 0] = beta(tab665_m[75, 0], tab665_m[146, 0])
tabk665[5, 1] = z(tab665_m[75, 0], tab665_m[146, 0])
tabk665[6, 0] = beta(tab665_m[146, 0], tab665_m[214, 0])
tabk665[6, 1] = z(tab665_m[146, 0], tab665_m[214, 0])
tabk665[7, 0] = beta(tab665_m[214, 0], tab665_m[275, 0])
tabk665[7, 1] = z(tab665_m[214, 0], tab665_m[275, 0])
tabk665[8, 0] = beta(tab665_m[275, 0], tab665_m[314, 0])
tabk665[8, 1] = z(tab665_m[275, 0], tab665_m[314, 0])

#Pi = 0.62 & Pr = 0 
plt.title('766')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab766_p[:, 0], tab766_p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab766_m[:, 0], tab766_m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab766_g[:, 0], tab766_g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

d766 = dist(13.0, 13.3)
id766 = idist(13.0, 13.3)

tabk766 = np.zeros((9,2))
tabk766[0, 0] = beta(tab766_p[43, 0],tab766_p[113, 0])
tabk766[0, 1] = z(tab766_p[43, 0], tab766_p[113, 0]) 
tabk766[1, 0] = beta(tab766_p[113, 0], tab766_p[187, 0])
tabk766[1, 1] = z(tab766_p[113, 0], tab766_p[187, 0])
tabk766[2, 0] = beta(tab766_p[187, 0], tab766_p[250, 0])
tabk766[2, 1] = z(tab766_p[187, 0], tab766_p[250, 0])
tabk766[3, 0] = beta(tab766_p[250, 0], tab766_p[297, 0])
tabk766[3, 1] = z(tab766_p[250, 0], tab766_p[297, 0])
tabk766[4, 0] = beta(tab766_m[16, 0], tab766_m[81, 0])
tabk766[4, 1] = z(tab766_m[16, 0], tab766_m[81, 0])
tabk766[5, 0] = beta(tab766_m[81, 0], tab766_m[152, 0])
tabk766[5, 1] = z(tab766_m[81, 0], tab766_m[152, 0])
tabk766[6, 0] = beta(tab766_m[152, 0], tab766_m[219, 0])
tabk766[6, 1] = z(tab766_m[152, 0], tab766_m[219, 0])
tabk766[7, 0] = beta(tab766_m[219, 0], tab766_m[282, 0])
tabk766[7, 1] = z(tab766_m[219, 0], tab766_m[282, 0])
tabk766[8, 0] = beta(tab766_m[282, 0], tab766_m[323, 0])
tabk766[8, 1] = z(tab766_m[282, 0], tab766_m[323, 0])

#Pi = 0.70 & Pr = 0.3
plt.title('880')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab880_p[:, 0], tab880_p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab880_m[:, 0], tab880_m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab880_g[:, 0], tab880_g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

d880 = dist(12.1, 12.2)
id880 = idist(12.1, 12.2)

tabk880 = np.zeros((9,2))
tabk880[0, 0] = beta(tab880_p[36, 0],tab880_p[105, 0])
tabk880[0, 1] = z(tab880_p[36, 0], tab880_p[105, 0]) 
tabk880[1, 0] = beta(tab880_p[105, 0], tab880_p[185, 0])
tabk880[1, 1] = z(tab880_p[105, 0], tab880_p[185, 0])
tabk880[2, 0] = beta(tab880_p[185, 0], tab880_p[247, 0])
tabk880[2, 1] = z(tab880_p[185, 0], tab880_p[247, 0])
tabk880[3, 0] = beta(tab880_p[247, 0], tab880_p[298, 0])
tabk880[3, 1] = z(tab880_p[247, 0], tab880_p[298, 0])
tabk880[4, 0] = beta(tab880_m[12, 0], tab880_m[75, 0])
tabk880[4, 1] = z(tab880_m[12, 0], tab880_m[75, 0])
tabk880[5, 0] = beta(tab880_m[75, 0], tab880_m[149, 0])
tabk880[5, 1] = z(tab880_m[75, 0], tab880_m[149, 0])
tabk880[6, 0] = beta(tab880_m[149, 0], tab880_m[219, 0])
tabk880[6, 1] = z(tab880_m[149, 0], tab880_m[219, 0])
tabk880[7, 0] = beta(tab880_m[219, 0], tab880_m[283, 0])
tabk880[7, 1] = z(tab880_m[219, 0], tab880_m[283, 0])
tabk880[8, 0] = beta(tab880_m[283, 0], tab880_m[324, 0])
tabk880[8, 1] = z(tab880_m[283, 0], tab880_m[324, 0])

#Pi = 0.66 & Pr = 0.2 A FAIRE
plt.title('1024')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab1024_p[:, 0], tab1024_p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab1024_m[:, 0], tab1024_m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab1024_g[:, 0], tab1024_g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

d1024 = dist(13.0, 14.0)
id1024 = idist(13.0, 14.0)

tabk1024 = np.zeros((9,2))
tabk1024[0, 0] = beta(tab1024_p[111, 0],tab1024_p[187, 0])
tabk1024[0, 1] = z(tab1024_p[111, 0], tab1024_p[187, 0]) 
tabk1024[1, 0] = beta(tab1024_p[187, 0], tab1024_p[251, 0])
tabk1024[1, 1] = z(tab1024_p[187, 0], tab1024_p[251, 0])
tabk1024[2, 0] = beta(tab1024_p[251, 0], tab1024_p[304, 0])
tabk1024[2, 1] = z(tab1024_p[251, 0], tab1024_p[304, 0])
tabk1024[3, 0] = beta(tab1024_p[304, 0], tab1024_p[334, 0])
tabk1024[3, 1] = z(tab1024_p[304, 0], tab1024_p[334, 0])
tabk1024[4, 0] = beta(tab1024_m[18, 0], tab1024_m[78, 0])
tabk1024[4, 1] = z(tab1024_m[18, 0], tab1024_m[78, 0])
tabk1024[5, 0] = beta(tab1024_m[78, 0], tab1024_m[153, 0])
tabk1024[5, 1] = z(tab1024_m[78, 0], tab1024_m[153, 0])
tabk1024[6, 0] = beta(tab1024_m[153, 0], tab1024_m[222, 0])
tabk1024[6, 1] = z(tab1024_m[153, 0], tab1024_m[222, 0])
tabk1024[7, 0] = beta(tab1024_m[222, 0], tab1024_m[282, 0])
tabk1024[7, 1] = z(tab1024_m[222, 0], tab1024_m[282, 0])
tabk1024[8, 0] = beta(tab1024_m[282, 0], tab1024_m[326, 0])
tabk1024[8, 1] = z(tab1024_m[282, 0], tab1024_m[326, 0])

#calcul des densités électroniques
