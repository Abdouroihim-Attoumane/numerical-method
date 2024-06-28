import numpy as np
import math
import matplotlib.pyplot as plt
import time
np.set_printoptions(precision=5) #n'affiche que 5 décimales pour les nombres décimaux

def Trace_Rect(F,a,b,N):
    Np=1000 #nombre de points placés sur chaque courbe.
    U=np.linspace(a,b,Np) #les U[j] sont les abscisses de tous les points que l'on va placer sur le graphique
    V=F(U); #les V[j] sont les ordonnées des points de la courbe de f
    #on stocke dans les tableaux ci-dessous les valeurs au points U[j] des fonctions par morceaux
    WG=np.zeros(Np)
    WPM=np.zeros(Np)
    WD=np.zeros(Np)
    X=np.linspace(a,b,N+1) #subdivision uniforme de [a,b]
    for j in range(Np):
        for i in range(N):
            if (X[i]<=U[j]) and (X[i+1]>U[j]):  #si U[j] appartient à l'intervalle [X[j],X[j+1][
                WG[j]=F(X[i])
                WD[j]=F(X[i+1])
                WPM[j]=F((X[i]+X[i+1])/2)
        #Le test suivant permet de compléter la dernier valeur du tableau W si nécessaire 
        if U[j]==X[N]:
            WG[j]=F(X[N-1])
            WD[j]=F(X[N])
            WPM[j]=F((X[N-1]+X[N])/2)
            
                
    #on trace la fonction f et les fonctions en escalier:
    plt.plot(U,V,label="Fonction F")
    plt.plot(U,WG,label="Rectangles à gauche")
    plt.plot(U,WD,label="Rectangles à droite")
    plt.plot(U,WPM,label="Rectangles point milieu")
    plt.legend()
    return []

def RG(f,a,b,N):
    x=np.linspace(a,b,N+1)
    s=0
    for i in range(N):
        s=s+f(x[i])
    return s*(b-a)/N

def RD(f,a,b,N):
    x=np.linspace(a,b,N+1)
    s=0
    for i in range(N):
        s=s+f(x[i+1])
    return s*(b-a)/N

def RPM(f,a,b,N):
    x=np.linspace(a,b,N+1)
    s=0
    for i in range(N):
        s=s+f((x[i]+x[i+1])/2)
    return s*(b-a)/N

def Trap(f,a,b,N):
    x=np.linspace(a,b,N+1)
    s=0
    for i in range(N):
        s=s+(f(x[i])+f(x[i+1]))/2
    return s*(b-a)/N

def Simpson(f,a,b,N):
    x=np.linspace(a,b,N+1)
    s=0
    for i in range(N):
        s=s+f(x[i])+4*f((x[i]+x[i+1])/2) + f(x[i+1])
    return s*(b-a)/(6*N)

def erreur(f,a,b,N,Iexact):
    #calcule pour chacune des méthodes des rectangles, l'erreur commise
    #entre l'intégrale calculée et la valeur exacte Iexact.
    #Les valeurs sont stockées dans un tableau E.
    E=np.zeros(5);
    E[0]=abs(RD(f,a,b,N)-Iexact)
    E[1]=abs(RG(f,a,b,N)-Iexact)
    E[2]=abs(RPM(f,a,b,N)-Iexact)
    E[3]=abs(Trap(f,a,b,N)-Iexact)
    E[4]=abs(Simpson(f,a,b,N)-Iexact)
    return E