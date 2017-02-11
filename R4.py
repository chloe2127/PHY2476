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
	return 1.113 - (d1+d2)/200

def idist(d1, d2): 				#incertitude sur la position (en m)
	return np.abs((d1-d2)/200)

def beta(l1, l2): 				#retourne k (beta) le vecteur d'onde (en m-1)
	l = np.abs(l1-l2)*2/100
	return 2*np.pi/l

def z(z1, z2): 					#position antenne (en cm)
	return (z1+z2)/200

def densite(tabK):
	# K : le tableau des k pour une certaine puissance
	# tabPhase : la colonne des k dans le tableau du digramma de phase

	densite = np.zeros((len(tabK[:,0]), 2)) # ce tableau va contenir toutes les densites calculees et l'incertitude associée

	for i in range(0,len(tabK[:,0])):
		idx = (np.abs(diagphase[:, 1] - tabK[i,0])).argmin()
		rapport = diagphase[idx, 1]
		densite[i, 0] = (omega**2)*me*epsilon0/((rapport**2)*(e**2))
		densite[i, 1] = np.abs((diagphase[idx+1, 1]-diagphase[idx, 1])/2) 	#incert (écart-type)

	return densite

def regression(x, y, exp, color): # Régression linéaire avec les formules 6.11,  6.23 du livre de Bevington
	
	# a,b = les paramètres du fit y=a+b*x; da, db, les incertitudes dans a et b
	N = x.size
	delta = N*np.sum(x**2) - (np.sum(x))**2
	a = (np.sum(x**2)*np.sum(y) - np.sum(x)*np.sum(x*y))/delta
	b = (N*np.sum(x*y) - np.sum(x)*np.sum(y))/delta
	s2 = np.sum((y-a-b*x)**2)/(N-2)
	da = np.sqrt(s2/delta*np.sum(x**2))
	db = np.sqrt(N*s2/delta)
	err = np.sqrt(s2)
	yf = a + b*x
	print(u'décalage=%.2f pente=%.2f \nsigma(décal.)=%.2f sigma(pente)=%.2f sigma=%.3f' % (a,b,da,db,err))
	
	plt.title(r'%s' %exp)
	plt.xlabel('z (m)')
	plt.ylabel(r'densité électronique ($m^{-3}$)')
	plt.plot(x, yf, color = r'%s' %color)
	#label = 'fit de %.3f' %err
	return err, b, a, da, db

def incert_theta(iPi, iPr, iz, ia, ib):
	return np.sqrt(iPi**2 + iPr**2 + iz**2 + ia**2 + ib**2)



###########################################################
# Définition des constantes

#constantes pour fct densité
omega = 2*np.pi*600e6			#2pi*f (600 MHz)
me = 9.10938356e-31 			#masse au repos d'un électron (kg)
epsilon0 = 8.8518782e-12		#permittivité du vide (m-3kg-1s4A2)
e = 1.60217662e-19*1e6			#charge électronique (C)

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
#plt.title("diagramme de phase de l'argon")
plt.xlabel(r'$\beta(cm^{-1})$')
plt.ylabel(r'$\omega/\omega_p$')
plt.scatter(diagphase[:, 1], diagphase[:, 0])
plt.show()

#1e5 sur les puissances (en mW)

#A PRESSION CONSTANTE (50 mTorrs)

#Pi = 0.21 & Pr = 0.02
#plt.title('1')
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
tabk1[:, 1] = d1 - tabk1[:, 1]

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
tabk2[:, 1] = d2 - tabk2[:, 1]

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
tabk3[:, 1] = d3 - tabk3[:, 1]

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
tabk4[:, 1] = d4 - tabk4[:, 1]

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
tabk5[:, 1] = d5 - tabk5[:, 1]

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
tabk6[:, 1] = d6 - tabk6[:, 1]

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
tabk7[:, 1] = d7 - tabk7[:, 1]

#calcul des densités électroniques

#plt.title('densité électronique')
plt.xlabel('position (m)')
plt.ylabel(r'densité électronique (m$^{-3}$)')
e1, b1, a1, siga1, sigb1 = regression(tabk1[:, 1], densite(tabk1)[:,0], '1', 'c')
plt.errorbar(tabk1[:,1], densite(tabk1)[:,0], xerr = id1, yerr = densite(tabk1)[:,1], fmt = '.', color = 'c', label='1')
e2, b2, a2, siga2, sigb2 = regression(tabk2[:, 1], densite(tabk2)[:,0], '2', 'b')
plt.errorbar(tabk2[:,1], densite(tabk2)[:,0], xerr = id2, yerr = densite(tabk2)[:,1], color = 'b', fmt = '.', label='2')
e3, b3, a3, siga3, sigb3 = regression(tabk3[:, 1], densite(tabk3)[:,0], '3', 'r')
plt.errorbar(tabk3[:,1], densite(tabk3)[:,0], xerr = id3, yerr = densite(tabk3)[:,1], color = 'r', fmt = '.', label='3')
e4, b4, a4, siga4, sigb4 = regression(tabk4[:, 1], densite(tabk4)[:,0], '4', 'g')
plt.errorbar(tabk4[:,1], densite(tabk4)[:,0], xerr = id4, yerr = densite(tabk4)[:,1], color = 'g', fmt = '.', label='4')
e5, b5, a5, siga5, sigb5 = regression(tabk5[:, 1], densite(tabk5)[:,0], '5', 'm')
plt.errorbar(tabk5[:,1], densite(tabk5)[:,0], xerr = id5, yerr = densite(tabk5)[:,1], color = 'm', fmt = '.', label='5')
e6, b6, a6, siga6, sigb6 = regression(tabk6[:, 1], densite(tabk6)[:,0], '1', 'y')
plt.errorbar(tabk6[:,1], densite(tabk6)[:,0], xerr = id6, yerr = densite(tabk6)[:,1], color = 'y', fmt = '.', label='6')
e7, b7, a7, siga7, sigb7 = regression(tabk7[:, 1], densite(tabk7)[:,0], '7', 'k')
plt.errorbar(tabk7[:,1], densite(tabk7)[:,0], xerr = id7, yerr = densite(tabk7)[:,1], color = 'k', fmt = '.', label='7')
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
tabk163[:, 1] = d163 - tabk163[:, 1]

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
tabk261[:, 1] = d261 - tabk261[:, 1]

#Pi = 0.72 & Pr = 0.01
plt.title('363')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab363_p[:, 0], tab363_p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab363_m[:, 0], tab363_m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab363_g[:, 0], tab363_g[:, 1], color = 'b', marker = '.', label='g')
plt.legend(loc='best')
plt.show()

d363 = dist(24.6, 24.2)
id363 = idist(24.6, 24.2)

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
tabk363[:, 1] = d363 - tabk363[:, 1]

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
tabk455[:, 1] = d455 - tabk455[:, 1]

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
tabk573[:, 1] = d573 - tabk573[:, 1]

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
tabk665[:, 1] = d665 - tabk665[:, 1]

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
tabk766[:, 1] = d766 - tabk766[:, 1]

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
tabk880[:, 1] = d880 - tabk880[:, 1]

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
tabk1024[:, 1] = d1024 - tabk1024[:, 1]

#calcul des densités électroniques
#plt.title('densité électronique')
plt.xlabel('position (m)')
plt.ylabel(r'densité électronique (m$^{-3}$)')
e163, b163, a163, siga163, sigb163 = regression(tabk163[:, 1], densite(tabk163)[:,0], '163', 'c')
plt.errorbar(tabk163[:, 1], densite(tabk163)[:,0], xerr = id163, yerr = densite(tabk163)[:,1], color = 'c', fmt = '.', label='163')
e261, b261, a261, siga261, sigb261 = regression(tabk261[:, 1], densite(tabk261)[:,0], '261', 'b')
plt.errorbar(tabk261[:,1], densite(tabk261)[:,0], xerr = id261, yerr = densite(tabk261)[:,1], color = 'b', fmt = '.', label='261')
e363, b363, a363, siga363, sigb363 = regression(tabk363[:, 1], densite(tabk363)[:,0], '363', 'r')
plt.errorbar(tabk363[:,1], densite(tabk363)[:,0], xerr = id363, yerr = densite(tabk363)[:,1], color = 'r', fmt = '.', label='363')
e455, b455, a455, siga455, sigb455 = regression(tabk455[:, 1], densite(tabk455)[:,0], '455', 'g')
plt.errorbar(tabk455[:, 1], densite(tabk455)[:,0], xerr = id455, yerr = densite(tabk455)[:,1], color = 'g', fmt = '.', label='455')
e573, b573, a573, siga573, sigb573 = regression(tabk573[:, 1], densite(tabk573)[:,0], '573', 'm')
plt.errorbar(tabk573[:,1], densite(tabk573)[:,0], xerr = id573, yerr = densite(tabk573)[:,1], color = 'm', fmt = '.', label='573')
e665, b665, a665, siga665, sigb665 = regression(tabk665[:, 1], densite(tabk665)[:,0], '665', 'y')
plt.errorbar(tabk665[:,1], densite(tabk665)[:,0], xerr = id665, yerr = densite(tabk665)[:,1], color = 'y', fmt = '.', label='665')
e766, b766, a766, siga766, sigb766 = regression(tabk766[:, 1], densite(tabk766)[:,0], '766', 'k')
plt.errorbar(tabk766[:,1], densite(tabk766)[:,0], xerr = id766, yerr = densite(tabk766)[:,1], color = 'k', fmt = '.', label='766')
e880, b880, a880, siga880, sigb880 = regression(tabk880[:, 1], densite(tabk880)[:,0], '880', 'k')
plt.errorbar(tabk880[:,1], densite(tabk880)[:,0], xerr = id880, yerr = densite(tabk880)[:,1], color = 'k', fmt = 'd', label='880')
e1024, b1024, a1024, siga1024, sigb1024 = regression(tabk1024[:, 1], densite(tabk1024)[:,0], '1024', 'k')
plt.errorbar(tabk1024[:,1], densite(tabk1024)[:,0], xerr = id1024, yerr = densite(tabk1024)[:,1], color = 'k', fmt = 'D', label='1024')
plt.legend(loc='best')
plt.show()

###################################################################################################
# longueur du plasma en fct de la puissance incidente
tabd = np.zeros((16, 2))
tabd[0,:] = d1, id1
tabd[1,:] = d2, id2
tabd[2,:] = d3, id3
tabd[3,:] = d4, id4
tabd[4,:] = d5, id5
tabd[5,:] = d6, id6
tabd[6,:] = d7, id7
tabd[7,:] = d163, id163
tabd[8,:] = d261, id261
tabd[9,:] = d363, id363
tabd[10,:] = d455, id455
tabd[11,:] = d573, id573
tabd[12,:] = d665, id665
tabd[13,:] = d766, id766
tabd[14,:] = d880, id880
tabd[15,:] = d1024, id1024

tabPuissInc = np.zeros((16, 2))
tabPuissInc[0,:] = .2e5, .01e5
tabPuissInc[1,:] = .245e5, .01e5
tabPuissInc[2,:] = .29e5, .01e5
tabPuissInc[3,:] = .35e5, .01e5
tabPuissInc[4,:] = .4e5, .01e5
tabPuissInc[5,:] = .47e5, .01e5
tabPuissInc[6,:] = .54e5, .01e5
tabPuissInc[7,:] = .56e5, .01e5
tabPuissInc[8,:] = .73e5, .01e5
tabPuissInc[9,:] = .72e5, .01e5
tabPuissInc[10,:] = .69e5, .01e5
tabPuissInc[11,:] = .7e5, .01e5
tabPuissInc[12,:] = .72e5, .01e5
tabPuissInc[13,:] = .62e5, .05e5
tabPuissInc[14,:] = .7e5, .01e5
tabPuissInc[15,:] = .66e5, .02e5

#plt.title('longueur en fct de Puissance (inc)')
plt.xlabel('Puissance incidente (mW)')
plt.ylabel(r"longeur du plasma L (m)")
plt.errorbar(tabPuissInc[:6, 0], tabd[:6, 0], xerr = tabPuissInc[:6, 1], yerr = tabd[:6, 1], fmt = '.')
plt.errorbar(tabPuissInc[7:, 0], tabd[7:, 0], xerr = tabPuissInc[7:, 1], yerr = tabd[7:, 1], fmt = '.')
plt.show()

####################################################################################################
# theta

tabPuissRef = np.zeros((16, 2))
tabPuissRef[0,:] = .02e5, .01e5
tabPuissRef[1,:] = .02e5, .01e5
tabPuissRef[2,:] = .025e5, .01e5
tabPuissRef[3,:] = .021e5, .01e5
tabPuissRef[4,:] = .026e5, .01e5
tabPuissRef[5,:] = .025e5, .01e5
tabPuissRef[6,:] = .028e5, .01e5
tabPuissRef[7,:] = .04e5, .01e5
tabPuissRef[8,:] = .02e5, .01e5
tabPuissRef[9,:] = .01e5, .01e5
tabPuissRef[10,:] = .01e5, .01e5
tabPuissRef[11,:] = .01e5, .01e5
tabPuissRef[12,:] = .04e5, .01e5
tabPuissRef[13,:] = .01e5, .01e5
tabPuissRef[14,:] = .03e5, .01e5
tabPuissRef[15,:] = .01e5, .01e5


# theta = Puissance absolue/(aL**2/2 - bL) (qui est le résultat de l'intégrale pour Ne)
tabTheta = np.zeros((16, 2))
tabTheta[0,:] = (tabPuissInc[0,0]-tabPuissRef[0,0])/(.5*a1*tabd[0,0]**2 + b1*tabd[0,0]), incert_theta(tabPuissInc[0, 1], tabPuissRef[0, 1], tabd[0, 1], siga1, sigb1)
tabTheta[1,:] = (tabPuissInc[1,0]-tabPuissRef[1,0])/(.5*a2*tabd[1,0]**2 + b2*tabd[1,0]), incert_theta(tabPuissInc[1, 1], tabPuissRef[1, 1], tabd[1, 1], siga2, sigb2)
tabTheta[2,:] = (tabPuissInc[2,0]-tabPuissRef[2,0])/(.5*a3*tabd[2,0]**2 + b3*tabd[2,0]), incert_theta(tabPuissInc[2, 1], tabPuissRef[2, 1], tabd[2, 1], siga3, sigb3)
tabTheta[3,:] = (tabPuissInc[3,0]-tabPuissRef[3,0])/(.5*a4*tabd[3,0]**2 + b4*tabd[3,0]), incert_theta(tabPuissInc[3, 1], tabPuissRef[3, 1], tabd[3, 1], siga4, sigb4)
tabTheta[4,:] = (tabPuissInc[4,0]-tabPuissRef[4,0])/(.5*a5*tabd[4,0]**2 + b5*tabd[4,0]), incert_theta(tabPuissInc[4, 1], tabPuissRef[4, 1], tabd[4, 1], siga5, sigb5)
tabTheta[5,:] = (tabPuissInc[5,0]-tabPuissRef[5,0])/(.5*a6*tabd[5,0]**2 + b6*tabd[5,0]), incert_theta(tabPuissInc[5, 1], tabPuissRef[5, 1], tabd[5, 1], siga6, sigb6)
tabTheta[6,:] = (tabPuissInc[6,0]-tabPuissRef[6,0])/(.5*a7*tabd[6,0]**2 + b7*tabd[6,0]), incert_theta(tabPuissInc[6, 1], tabPuissRef[6, 1], tabd[6, 1], siga7, sigb7)
tabTheta[7,:] = (tabPuissInc[7,0]-tabPuissRef[7,0])/(.5*a163*tabd[7,0]**2 + b163*tabd[7,0]), incert_theta(tabPuissInc[7, 1], tabPuissRef[7, 1], tabd[7, 1], siga163, sigb163)
tabTheta[8,:] = (tabPuissInc[8,0]-tabPuissRef[8,0])/(.5*a261*tabd[8,0]**2 + b261*tabd[8,0]), incert_theta(tabPuissInc[8, 1], tabPuissRef[8, 1], tabd[8, 1], siga261, sigb261)
tabTheta[9,:] = (tabPuissInc[9,0]-tabPuissRef[9,0])/(.5*a363*tabd[9,0]**2 + b363*tabd[9,0]), incert_theta(tabPuissInc[9, 1], tabPuissRef[9, 1], tabd[9, 1], siga363, sigb363)
tabTheta[10,:] = (tabPuissInc[10,0]-tabPuissRef[10,0])/(.5*a455*tabd[10,0]**2 + b455*tabd[10,0]), incert_theta(tabPuissInc[10, 1], tabPuissRef[10, 1], tabd[10, 1], siga455, sigb455)
tabTheta[11,:] = (tabPuissInc[11,0]-tabPuissRef[11,0])/(.5*a573*tabd[11,0]**2 + b573*tabd[11,0]), incert_theta(tabPuissInc[11, 1], tabPuissRef[11, 1], tabd[11, 1], siga573, sigb573)
tabTheta[12,:] = (tabPuissInc[12,0]-tabPuissRef[12,0])/(.5*a665*tabd[12,0]**2 + b665*tabd[12,0]), incert_theta(tabPuissInc[12, 1], tabPuissRef[12, 1], tabd[12, 1], siga665, sigb665)
tabTheta[13,:] = (tabPuissInc[13,0]-tabPuissRef[13,0])/(.5*a766*tabd[13,0]**2 + b766*tabd[13,0]), incert_theta(tabPuissInc[13, 1], tabPuissRef[13, 1], tabd[13, 1], siga766, sigb766)
tabTheta[14,:] = (tabPuissInc[14,0]-tabPuissRef[14,0])/(.5*a880*tabd[14,0]**2 + b880*tabd[14,0]), incert_theta(tabPuissInc[14, 1], tabPuissRef[14, 1], tabd[14, 1], siga880, sigb880)
tabTheta[15,:] = (tabPuissInc[15,0]-tabPuissRef[15,0])/(.5*a1024*tabd[15,0]**2 + b1024*tabd[15,0]), incert_theta(tabPuissInc[15, 1], tabPuissRef[15, 1], tabd[15, 1], siga1024, sigb1024)
print(tabTheta)

