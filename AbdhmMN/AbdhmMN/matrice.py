#Librairies nécessaires:

import numpy as np

# La matrice de Vandermonde sans boucle
#x doit être de taille (n,0): x=[a,b,c, ...] dans les quatre fonctions ci-dessous.

def vand(x):
    n=np.size(x) #taille de la matrice A
    #on convertit x en une matrice X de taille (1,n)
    X=np.zeros((1,n))
    X[0,:]=x 
    #Construction de la matrice A
    A=np.dot(np.ones((n,1)),X)
    #A est la matrice n-n dont chaque ligne est le vecteur x
    #Il faut que X soit un vecteur ligne (1,n)!
    z=np.zeros((1,n)) #on impose que z est un vecteur ligne
    z[:]=np.arange(n) #z=[0 1 2 ... n-1]
    B=np.transpose(np.dot(np.ones((n,1)),z))
    #B est la matrice dont chaque colonne est le vecteur [0 1 ... (n-1)]
    A=A**B
    #cette commande permet de caculer les puissances terme à terme
    return A
