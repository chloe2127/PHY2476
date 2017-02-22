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

def comparaison(noir, tab, Pmax, dim, Id):			#paramètres à comparer pour grosse vs petite
	Is = noir[:,0]/(np.exp((q*noir[:,1])/k*T)-1)	#courant de saturation qu'on calcule dans le noir (pour I=Id)
	phi = Pmax/(dim[0,0]*dim[1,0])					#phi = Pmax/surface
	plt.scatter(tab[:,1], np.log(Id[:]))
	plt.show()										#juste pour vérifier que le choix des pts pour la pente et correct
	n = (tab[-1,1]-tab[-10,1])/(np.log(Id[-1])-np.log(Id[-10])) #inverse pente ln(Id) vs V

	return Is, phi, n

def magie(tab):
	
	iIcc = tab[:,1].searchsorted(0)			 	#trouve l'indice où I(V=0)
	iIcc2 = tab[:,1].searchsorted(tab[iIcc,1])	#deuxième indice proche de V=0
	Icc = (tab[iIcc,0]+tab[iIcc2,0])/2
	
	Id = tab[:,0] - Icc

	iVco = tab[:,0].searchsorted(0)			 	#trouve l'indice où V(I=0)
	iVco2 = tab[:,0].searchsorted(tab[iVco,0])	#deuxième indice proche de I=0
	Vco = (tab[iVco,1]+tab[iVco2,1])/2
	
	Rp = (tab[iIcc+1,1]-tab[iIcc-1,1])/(tab[iIcc+1,0]-tab[iIcc-1,0])	#inverse pente tangente en Icc
	Rs = (tab[iVco+1,1]-tab[iVco-1,1])/(tab[iVco+1,0]-tab[iVco-1,0])	#inverse pente tangente en Icc

	tabPmax = np.zeros(len(tab[:,1]))			#array où on va mettre les Pmax provisoirs
	for i in range(len(tab[:,0])):				#trouver Pmax
		tabPmax[i] = np.amax(np.abs(tab[i,0])*np.abs(tab[:,1])) #pour ce I, Pmax (pour tous les V)
		
	iPmax = np.argmax(tabPmax) 					#indice de Pmax pour I et V confondus
	Pmax = tabPmax[iPmax]						
	Ipm = tab[iPmax,0]							#I qui permet Pmax
	Vpm = Pmax/Ipm 								#V qui permet Pmax

	vals = np.array([Icc, Vco, Rp, Rs, Pmax, Ipm, Vpm])
	
	return vals, Id


###########################################################
# Définition des constantes

k = 1.38e-23	#cst de boltzmann [J/K]
T = 300			#température de la cellule [K] = temperature ambiante (on l'utilise pour calculer Is dans le noir)
q = 1.6e-19		#charge élémentaire [C] q=e/n

dimp = np.array([[(3.9+4.2)/2, (4.2-3.9)/2], [(1.9+1.5)/2, (1.9-1.5)/2]]) 		#dimensions et incertitude petite
dimg = np.array([[(11.1+10.9)/2, (11.1-10.9)/2], [(6.9+6.4)/2, (6.9-6.4)/2]])	#dimensions et incertitude grande

#LECTURE DES DONNÉES
#petites cellules (en série et en parallèle)
pnoir = np.loadtxt('black.csv', delimiter=';', usecols=range(2))

tabps1_1 = np.loadtxt('ps1_1.csv', delimiter=';', usecols=range(2))
tabps1_2 = np.loadtxt('ps1_2.csv', delimiter=';', usecols=range(2))
tabps1_3 = np.loadtxt('ps1_3.csv', delimiter=';', usecols=range(2))
tabps1_4 = np.loadtxt('ps1_4.csv', delimiter=';', usecols=range(2))
#tabps1 = np.concatenate((tabps1_1, tabps1_2, tabps1_3, tabps1_4), axis=1)
#print(tabps1)

tabps2_1 = np.loadtxt('ps2_1.csv', delimiter=';', usecols=range(2))
tabps2_2 = np.loadtxt('ps2_2.csv', delimiter=';', usecols=range(2))
tabps2_3 = np.loadtxt('ps2_3.csv', delimiter=';', usecols=range(2))
tabps2_4 = np.loadtxt('ps2_4.csv', delimiter=';', usecols=range(2))

tabps3_1 = np.loadtxt('ps3_1.csv', delimiter=';', usecols=range(2))
tabps3_2 = np.loadtxt('ps3_2.csv', delimiter=';', usecols=range(2))
tabps3_3 = np.loadtxt('ps3_3.csv', delimiter=';', usecols=range(2))
tabps3_4 = np.loadtxt('ps3_4.csv', delimiter=';', usecols=range(2))

tabps4_1 = np.loadtxt('ps4_1.csv', delimiter=';', usecols=range(2))
tabps4_2 = np.loadtxt('ps4_2.csv', delimiter=';', usecols=range(2))
tabps4_3 = np.loadtxt('ps4_3.csv', delimiter=';', usecols=range(2))
tabps4_4 = np.loadtxt('ps4_4.csv', delimiter=';', usecols=range(2))

tabpp2_1 = np.loadtxt('pp2_1.csv', delimiter=';', usecols=range(2))
tabpp2_2 = np.loadtxt('pp2_2.csv', delimiter=';', usecols=range(2))
tabpp2_3 = np.loadtxt('pp2_3.csv', delimiter=';', usecols=range(2))
tabpp2_4 = np.loadtxt('pp2_4.csv', delimiter=';', usecols=range(2))

tabpp3_1 = np.loadtxt('pp3_1.csv', delimiter=';', usecols=range(2))
tabpp3_2 = np.loadtxt('pp3_2.csv', delimiter=';', usecols=range(2))
tabpp3_3 = np.loadtxt('pp3_3.csv', delimiter=';', usecols=range(2))
tabpp3_4 = np.loadtxt('pp3_4.csv', delimiter=';', usecols=range(2))

tabpp4_1 = np.loadtxt('pp4_1.csv', delimiter=';', usecols=range(2))
tabpp4_2 = np.loadtxt('pp4_2.csv', delimiter=';', usecols=range(2))
tabpp4_3 = np.loadtxt('pp4_3.csv', delimiter=';', usecols=range(2))
tabpp4_4 = np.loadtxt('pp4_4.csv', delimiter=';', usecols=range(2))

#grosse cellule
gnoir = np.loadtxt('gwo_black_bandit.csv', delimiter=';', usecols=range(2))

tabg1 = np.loadtxt('g_1.csv', delimiter=';', usecols=range(2))
print(len(tabg1))
tabg2 = np.loadtxt('g_2.csv', delimiter=';', usecols=range(2))
tabg3 = np.loadtxt('g_3.csv', delimiter=';', usecols=range(2))
tabg4 = np.loadtxt('g_4.csv', delimiter=';', usecols=range(2))

###########################################################
# Progamme principal

#PETITE(S)

#ps1:
valsps1_1, Idps1_1 = magie(tabps1_1)
print(valsps1_1)
valsps1_2, Idps1_2 = magie(tabps1_2)
valsps1_3, Idps1_3 = magie(tabps1_3)
valsps1_4, Idps1_4 = magie(tabps1_4)

#ps2:
valsps2_1, Idps2_1 = magie(tabps2_1)
valsps2_2, Idps2_2 = magie(tabps2_2)
valsps2_3, Idps2_3 = magie(tabps2_3)
valsps2_4, Idps2_4 = magie(tabps2_4)

#ps3:
valsps3_1, Idps3_1 = magie(tabps3_1)
valsps3_2, Idps3_2 = magie(tabps3_2)
valsps3_3, Idps3_3 = magie(tabps3_3)
valsps3_4, Idps3_4 = magie(tabps3_4)

#ps4:
valsps4_1, Idps4_1 = magie(tabps4_1)
valsps4_2, Idps4_2 = magie(tabps4_2)
valsps4_3, Idps4_3 = magie(tabps4_3)
valsps4_4, Idps4_4 = magie(tabps4_4)

#pp2:
valspp2_1, Idpp2_1 = magie(tabpp2_1)
valspp2_2, Idpp2_2 = magie(tabpp2_2)
valspp2_3, Idpp2_3 = magie(tabpp2_3)
valspp2_4, Idpp2_4 = magie(tabpp2_4)

#pp3:
valspp3_1, Idpp3_1 = magie(tabpp3_1)
valspp3_2, Idpp3_2 = magie(tabpp3_2)
valspp3_3, Idpp3_3 = magie(tabpp3_3)
valspp3_4, Idpp3_4 = magie(tabpp3_4)

#ps4:
valspp4_1, Idpp4_1 = magie(tabpp4_1)
valspp4_2, Idpp4_2 = magie(tabpp4_2)
valspp4_3, Idpp4_3 = magie(tabpp4_3)
valspp4_4, Idpp4_4 = magie(tabpp4_4)

#noir:
isp, phip, np = comparaison(pnoir, tabps1_1, valsps1_1[4], dimp, Idps1_1)

#GROSSE

#g1:
#valsg1, Idg1 = magie(tabg1)

#g2:
#valsg2, Idg2 = magie(tabg2)

#g3:
#valsg3, Idg3 = magie(tabg3)

#g4:
#valsg4, Idg4 = magie(tabg4)

#noir:
#isg, phig, ng = comparaison(gnoir, tabg1, valsg1[4], dimp, Idg1)



