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
# Définition des constantes

# lecture des données txt
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

#1e5 sur les puissances (en mW)

#Pi = 0.21 & Pr = 0.02
plt.title('1')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab1p[:, 0], tab1p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab1m[:, 0], tab1m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab1g[:, 0], tab1g[:, 1], color = 'b', marker = '.', label='g')
plt.show()

#Pi = 0.245 & Pr = 0.02
plt.title('2')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab2p[:, 0], tab2p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab2m[:, 0], tab2m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab2g[:, 0], tab2g[:, 1], color = 'b', marker = '.', label='g')
plt.show()

#Pi = 0.290 & Pr = 0.025
plt.title('3')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab3p[:, 0], tab3p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab3m[:, 0], tab3m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab3g[:, 0], tab3g[:, 1], color = 'b', marker = '.', label='g')
plt.show()

#Pi = 0.35 & Pr = 0.021
plt.title('4')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab4p[:, 0], tab4p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab4m[:, 0], tab4m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab4g[:, 0], tab4g[:, 1], color = 'b', marker = '.', label='g')
plt.show()

#Pi = 0.40 & Pr = 0.026
plt.title('5')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab5p[:, 0], tab5p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab5m[:, 0], tab5m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab5g[:, 0], tab5g[:, 1], color = 'b', marker = '.', label='g')
plt.show()

#Pi = 0.47 & Pr = 0.025
plt.title('6')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab6p[:, 0], tab6p[:, 1], color = 'c', marker = '.', label='p')
plt.scatter(tab6m[:, 0], tab6m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab6g[:, 0], tab6g[:, 1], color = 'b', marker = '.', label='g')
plt.show()

#Pi = 0.54 & Pr = 0.028
plt.title('7')
plt.xlabel('position (cm)')
plt.ylabel('amplitude (u.A.)')
plt.scatter(tab7m[:, 0], tab7m[:, 1], color = 'm', marker = '.', label='m')
plt.scatter(tab7g[:, 0], tab7g[:, 1], color = 'b', marker = '.', label='g')
plt.show()