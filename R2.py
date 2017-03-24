# -*- coding: utf-8 -*-

###########################################################
# PHY2476 - 2 Analyse de Fourier
# Chloé Lefebvre, Olivier Petit Vincent
# 21 mars 2017
###########################################################

###########################################################
# Importation des modules et des fonctions externes
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

###########################################################
# Définition de fonctions locales

def fourier(tab, forme):

	if forme == 'Sinus':
		plt.plot(t, np.sin(2*np.pi*f*t), color='k')
		#fctexp = np.empty(len(t))		#future fct d'onde
		fctexp = tab[0,1]*np.sqrt(2)*np.cos(tab[0,0]*omega*t + np.radians(tab[0,3]))
		print(fctexp)
		# for i in range(len(tab)): 		#boucle sur chaque n
		# 	if tab[i,3]<0: 				#si phin<0, phin = 0
		# 		fctexp += tab[i,1]*np.sqrt(2)*np.cos(tab[i,0]*omega*t)
		# 	else: 
		# 		fctexp += tab[i,1]*np.sqrt(2)*np.cos(tab[i,0]*omega*t + np.radians(tab[i,3]))

	elif forme == 'Carrée':
		plt.plot(t, (-1)*signal.square(2*np.pi*f*t), color='k')
		fctexp = np.empty(len(t))	
		for i in range(len(tab)): 		
			if (tab[i,0]%2) == 0: 	#on enleve les n pairs pour la fct carree
				fctexp += 0
				print('-pairs', fctexp)
			else: 
				print('phi', tab[i,3])
				if tab[i,3]<0: 				#si phin<0, phin = 0
					fctexp += tab[i,1]*np.sqrt(2)*np.cos(tab[i,0]*omega*t)
					print('-phin neg', fctexp)
				else: 
					fctexp += tab[i,1]*np.sqrt(2)*np.cos(tab[i,0]*omega*t + np.radians(tab[i,3]))
					print('normal', fctexp)
					

	elif forme == 'Triangulaire symétrique':
		plt.plot(t, (-1)*signal.sawtooth(2*np.pi*f*t, .5), color='k')	#0.5 for triangle wave, *(-1) pour la retourner (mettre par dessus nos donnees)
		fctexp = np.empty(len(t))		#future fct d'onde 
		for i in range(len(tab)): 		#boucle sur chaque n  
			if (tab[i,0]%2) == 0: 	#on enleve les n pairs pour la fct triangulaire
				fctexp += 0
			else: 
				if tab[i,3]<0: 				#si phin<0, phin = 0
					fctexp += tab[i,1]*np.sqrt(2)*np.cos(tab[i,0]*omega*t)
				else: 
					fctexp += tab[i,1]*np.sqrt(2)*np.cos(tab[i,0]*omega*t + np.radians(tab[i,3]))
			
	elif forme == 'Triangulaire asymétrique':
		plt.plot(t, signal.sawtooth(2*np.pi*f*t, 0), color='k')		#0 for falling ramp, 1 for rising ramp
		fctexp = np.empty(len(t))		#future fct d'onde 
		for i in range(len(tab)): 		#boucle sur chaque n  
			if (tab[i,0]%2) == 0: 	#on enleve les n pairs pour la fct triangulaire
				fctexp += 0
			else: 
				if tab[i,3]<=0: 				#si phin<0, phin = 0
					fctexp += tab[i,1]*np.sqrt(2)*np.cos(tab[i,0]*omega*t)
				else: 
					fctexp += tab[i,1]*np.sqrt(2)*np.cos(tab[i,0]*omega*t + np.radians(tab[i,3]))
			
	else: 						#forme == circuit RC
		fctexp = np.empty(len(t))		#future fct d'onde
		for i in range(len(tab)): 		#boucle sur chaque n
			if (tab[i,0]%2) == 0: 	#on enleve les n pairs pour la fct carree MEME ATTENUEE?
				fctexp += 0
			else: 
				if tab[i,3]<0: 				#si phin<0, phin = 0
					fctexp += tab[i,1]*np.sqrt(2)*np.cos(tab[i,0]*omega*t)
				else: 
					fctexp += tab[i,1]*np.sqrt(2)*np.cos(tab[i,0]*omega*t + np.radians(tab[i,3]))

	plt.scatter(t, fctexp/np.amax(tab[:,1]*np.sqrt(2)), label='%s' %forme)	#normaliser l'axe y pour voir lien avec la fct theo
	plt.xlabel('temps (s)')
	plt.xlim(0, 2/f)
	plt.ylabel('f(t)')
	plt.legend(loc='best')
	plt.show()

	return fctexp

def param(tab, forme):

	fig= plt.figure(1)
	
	plt.errorbar(tab[:,0], tab[:,1]*np.sqrt(2), yerr=tab[:,2]*np.sqrt(2), fmt='o', color='g')
	plt.ylabel(r'$C_n$ (mV)')
	plt.yticks(color='g')

	plt.twinx()
	plt.errorbar(tab[:,0], tab[:,3], yerr=tab[:,4], fmt='d', color='m', label=r'%s' %forme)
	plt.ylabel(r'$\phi_n (degrés)$')
	plt.ylim(-20, 200)
	plt.yticks(color='m')

	plt.xlabel('Harmonique')
	plt.xlim(tab[0,0]-1, tab[-1,0]+1)
	plt.legend(loc='upper right')
	plt.show()

	return fig

###########################################################
# Définition de constantes

f = 1e3				#f = 1 kHz
t = np.linspace(0, 2/f, 100, endpoint = True)	#2 periodes
omega = 2*np.pi*f	
R = 5e4				#ohm
C = 2.2e-6			#farad
forme = np.array(['Sinus', 'Carrée', 'Triangulaire asymétrique', 'Triangulaire symétrique', 'circuit RC'])

# lecture des fichiers -- [:,0] = harmonique n, [:,1] = cn/sqrt(2) (mV), [:,2] = incert cn (mV), [:,3] = phin (deg), [:,4] = incert phin (deg)
tabSin = np.loadtxt('sinus.txt') 
tabSquare = np.loadtxt('carre.txt') 
tabTriAsym = np.loadtxt('tri_asym.txt') 
tabTriSym = np.loadtxt('tri_sym.txt') 
tabRC = np.loadtxt('RC.txt')

###########################################################
# Programme principal

#graphiques de cn et phin
#figSin = param(tabSin, forme[0])
#figSquare = param(tabSquare, forme[1])
#figTriAsym = param(tabTriAsym, forme[2])
#figTriSym = param(tabTriSym, forme[3])
#figRC = param(tabRC, forme[4])

#graphique de transfo de fourier
fctSin = fourier(tabSin, forme[0])
fctSquare = fourier(tabSquare, forme[1])
fctTriAsym = fourier(tabTriAsym, forme[2])
fctTriSym = fourier(tabTriSym, forme[3])
#trouver la fct theo 
fctRC = fourier(tabRC, forme[4])
