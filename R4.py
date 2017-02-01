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
