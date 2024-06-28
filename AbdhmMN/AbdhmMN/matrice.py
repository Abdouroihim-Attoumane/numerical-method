#Librairies nécessaires:

import numpy as np

#--------------------------------------------------------------

#La fonction Tri-diagonal(n) construit la matrice diagonale demandée, de taille nxn

def TriD(taille, diag_valeur, sousdiag_valeur, surdiag_valeur):

    A=np.zeros((taille,taille))
    A[taille-1,taille-1] =  diag_valeur

    for i in range(n-1):
        A[i,i] = diag_valeur
        A[i+1,i] = sousdiag_valeur
        A[i,i+1] = surdiag_valeur

    return A

#--------------------------------------------------------------

#Permutation des colonnes de la matrice M les une à près les autres 

def PermutCol(M):

    (m,n)=np.shape(M);
    #on détermine le nombre de lignes et de colonnes de M
    N=np.zeros((m,n));#on crée une matrice de même taille
    #On remplit ensuite les colonnes de cette nouvelle matrice:
    N[0:m,0:n-1]=M[0:m,1:n]
    N[0:m,n-1]=M[0:m,0]
    return N

#--------------------------------------------------------------

# La matrice de Vandermonde sans boucle
#le vecteur doit être de taille (n,0): x=[a,b,c, ...] dans les quatre fonctions ci-dessous.

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
