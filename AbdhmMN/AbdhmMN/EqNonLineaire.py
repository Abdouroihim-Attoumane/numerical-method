#Librairies nécessaires

import numpy as np
import math
import matplotlib.pyplot as plt

#Fonctions Dicho, Secante et Newton  

def Dicho(f,a,b,N):
    A=np.zeros(N+1)
    B=np.zeros(N+1)
    X=np.zeros(N+1)
    if f(a)*f(b)>0:
        print("La fonction doit changer de signe entre a et b.")
        return []
    else:
        A[0]=a
        B[0]=b
        X[0]=(a+b)/2
        for k in range(N):
            if f(A[k])*f(X[k]) <= 0:
                A[k+1]=A[k]
                B[k+1]=X[k]
                X[k+1]=(A[k+1]+B[k+1])/2
            else:
                A[k+1]=X[k]
                B[k+1]=B[k]
                X[k+1]=(A[k+1]+B[k+1])/2
        return X

def Secante(f,a,b,N):
    X=np.zeros(N+1)
    X[0]=a;
    X[1]=b;
    for k in range(1,N):
        if X[k]==X[k-1]:
            X[k+1]=X[k]
        else:
            X[k+1]=X[k] - f(X[k])*(X[k]-X[k-1])/(f(X[k])-f(X[k-1]))
    return X

def Newton(f,g,a,b,N):
    X=np.zeros(N+1)
    X[0]=(a+b)/2
    for k in range(N):
        X[k+1]=X[k] - f(X[k])/g(X[k])
    return X

###################################
#Fonctions Nb_iterations

def NbIterations_Dicho(f,a,b,epsilon,s):
    N=2
    x=Dicho(f,a,b,N)
    while abs(x[N]-s)>epsilon:
        N=N+1
        x=Dicho(f,a,b,N)
    return N

def NbIterations_Secante(f,a,b,epsilon,s):
    N=2
    x=Secante(f,a,b,N)
    while abs(x[N]-s)>epsilon:
        N=N+1
        x=Secante(f,a,b,N)
    return N

def NbIterations_Newton(f,g,a,b,epsilon,s):
    N=2
    x=Newton(f,g,a,b,N)
    while abs(x[N]-s)>epsilon:
        N=N+1
        x=Newton(f,g,a,b,N)
    return N

###################################
#Fonctions Tracer

def Tracer_Dicho(f,a,b,N,s):
    E=Erreur_Dicho(f,a,b,N,s)
    plt.plot(-np.log(E[0:N]),-np.log(E[1:N+1]))
    plt.plot(-np.log(E[0:N]),-np.log(E[0:N])) #on place les points de même abscisse et sur la droite y=x
    plt.title("Evolution de e_(n+1) en fonction de e_n\n Méthode par dichotomie")
    plt.axis("equal")
    plt.xlabel("-log(e_n)")
    plt.ylabel("-log(e_(n+1))")
    plt.legend(("e_(n+1) en fonction de e_n","y=x"),loc="best")
    plt.grid() #pour afficher une grille
    plt.show()
    
def Tracer_Secante(f,a,b,N,s):
    E=Erreur_Secante(f,a,b,N,s)
    plt.plot(-np.log(E[0:N]),-np.log(E[1:N+1]))
    plt.title("Evolution de e_(n+1) en fonction de e_n\n Méthode de la sécante")
    plt.axis("equal")
    plt.xlabel("-log(e_n)")
    plt.ylabel("-log(e_(n+1))")
    plt.grid() #pour afficher une grille
    plt.show()

def Tracer_Newton(f,g,a,b,N,s):
    E=Erreur_Newton(f,g,a,b,N,s)
    plt.plot(-np.log(E[0:N]),-np.log(E[1:N+1]))
    plt.title("Evolution de e_(n+1) en fonction de e_n\n Méthode de Newton")
    plt.axis("equal")
    plt.xlabel("-log(e_n)")
    plt.ylabel("-log(e_(n+1))")
    plt.grid() #pour afficher une grille
    plt.show()