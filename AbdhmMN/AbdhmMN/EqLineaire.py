import numpy as np
import time
import math
import matplotlib.pyplot as plt

#Méthode de Gauss

def echangeligne(A,i,j):
    (n,p)=np.shape(A)
    E=np.zeros((n,p))
    E[:,:]=A[:,:]
    E[i-1,:]=A[j-1,:]
    E[j-1,:]=A[i-1,:]
    return E

def modifLj(A,j,i):  
    #Lj<-Lj - aji/aii Li
    (n,p)=np.shape(A)
    E=np.zeros((n,p))
    E[:,:]=A[:,:]
    E[j-1,:]=E[j-1,:]-(E[j-1,i-1]/E[i-1,i-1])*E[i-1,:]
    return E

def pivotmax(A,j):
    Lmax=j
    a=abs(A[j-1,j-1])
    (n,p)=np.shape(A)
    for k in range(j+1,n+1):
        if abs(A[k-1,j-1])>a:  
            Lmax=k
            a=abs(A[k-1,j-1])
    return Lmax

def TrigResolSup(T,C):
    (n,p)=np.shape(T)
    X=np.zeros((n,1))
    X[n-1]=C[n-1]/T[n-1,n-1]
    for j in range(n-1,0,-1):  #j va de n-1 à 0 avec un pas de -1
        s=0
        for k in range(j+1,n+1):
            s=s+T[j-1,k-1]*X[k-1]
            X[j-1] =(C[j-1]-s)/T[j-1,j-1]
    return X

def triangle(M):
    (n,p)=np.shape(M)
    EM=np.zeros((n,p))
    EM[:,:]=M[:,:]
    for i in range(1,n):
        k=pivotmax(EM,i)
        T=echangeligne(EM,i,k)
        for j in range(i+1,n+1): 
            EM=modifLj(EM,j,i)
    return EM

def Gauss(A,B): #B doit être un vecteur (n,1)
    (n,p)=np.shape(A)
    AB=np.zeros((n,p+1))
    AB[:,0:p]=A[:,0:p]
    AB[:,p]=B[:,0]
    EM=triangle(AB)
    T=EM[0:n,0:p]
    C=EM[0:n,p]
    X=TrigResolSup(T,C)
    return X

# Factorisation LU

def Transvection(n,i,j,a): 
#matrice de transvection Tij(a)
#Tij(a)*A revient à faire Li <- Li +aLj dans A
    T=np.zeros((n,n))
    for k in range(1,n+1):
        T[k-1,k-1]=1
    T[i-1,j-1]=a
    return T

def TrigResolInf(L,B): 
    #resolution de LY=B avec L triangulaire inférieure et inversible
    #on utilise la formule:
    #Y[j] = (B[j] -L[j,1] Y[1] - L[j,2] Y[2] - ... -L[j,j-1] Y[j-1])/L[j,j]
    (n,p)=np.shape(L)
    Y=np.zeros((n,1))
    Y[0]=B[0]/L[0,0]
    for j in range(2,n+1):
        s=0
        for k in range(1,j):
            s=s+L[j-1,k-1]*Y[k-1]
        Y[j-1]=(B[j-1]-s)/L[j-1,j-1]
    return Y
    
def InverseTrig(L):
    #calcule l'inverse de L avec L triangulaire inférieure
    #on note IL l'inverse de L. L'idée est IL=[F1, ... , Fn] avec IL*Ei=Fi où (E1,...,En) base canonique. 
    #Donc Fi est l'unique solution du système LX=Ei
    #on utilise la fonction triginf pour calculer Fi
    (n,p)=np.shape(L)
    IL=np.zeros((n,n))
    for i in range(1,n+1):
        E=np.zeros((n,1))
        E[i-1,0]=1
        IL[:,i-1]=TrigResolInf(L,E)[:,0]
    return IL

def LU(A): #décomposition LU de A
    (n,p)=np.shape(A)
    L=np.zeros((n,p))
    U=np.zeros((n,p))
    if n!=p: 
        print("Erreur, la matrice A n est pas carrée. \n")
    else:
        U[:,:]=A[:,:]
        IL=np.zeros((n,n))
        for j in range(1,n+1):
            IL[j-1,j-1]=1
        for i in range(1,n):
            for j in range(i+1,n+1):
                if U[i-1,i-1]==0:
                    print("Erreur, la matrice A ne possède pas de décomposition LU.")
                else:
                    a=U[j-1,i-1]/U[i-1,i-1]
                    T=Transvection(n,j,i,-a)
                    U=np.dot(T,U)
                    IL=np.dot(T,IL)
    L=InverseTrig(IL)
    return (L,U)

def ResolLU(A,B): 
    #résolution de AX=B avec la décomposition LU de A
    (L,U)=LU(A)
    Y=TrigResolInf(L,B)
    X=TrigResolSup(U,Y)
    return X

#Tests sur des exemples de grande taille

def TriDiag(n):
    A=np.zeros((n,n))
    for i in range(0,n-1):
        A[i,i]=3
        A[i,i+1]=-1
        A[i+1,i]=-2
    A[n-1,n-1]=3
    return A
    
def TracerTemps(N):
    tempsG=np.zeros(N)
    tempsLU=np.zeros(N)
    
    #Mesure des temps
    for n in range(1,N+1):
        A=TriDiag(n)
        XT=np.ones((n,1))
        B=np.dot(A,XT)
        #Gauss
        tic=time.time()
        Gauss(A,B)
        toc=time.time()
        tempsG[n-1]=toc-tic
        #LU
        tic=time.time()
        ResolLU(A,B)
        toc=time.time()
        tempsLU[n-1]=toc-tic
    
    #Courbes
    plt.plot(np.arange(1,N+1),tempsG,label="Gauss")
    plt.plot(np.arange(1,N+1),tempsLU,label="LU")
    plt.xlabel("Nombre d équation n")
    plt.ylabel("Temps t_n")
    plt.legend()
    plt.title("Evolution du temps de calcul en fonction du nombre d'équations")
    plt.show()
    return []

#Mesure de l'erreur commise

def Norme(X):
    #X est un vecteur colonne
    U=math.sqrt(np.dot(np.transpose(X),X)) #U est un vecteur (1,1)
    return U

def TracerErreur1(N):
    eG=np.zeros(N)
    eP=np.zeros(N)
    eL=np.zeros(N)
    
    #Erreur
    for n in range(1,N+1):
        A=TriDiag(n);#A est une matrice tridiagonale de rang n
        XT=np.ones((n,1))
        B=np.dot(A,XT)
        eG[n-1]=Norme(Gauss(A,B)-XT)
        eP[n-1]=Norme(np.linalg.solve(A,B)-XT)
        eL[n-1]=Norme(ResolLU(A,B)-XT)
    
    #Courbes
    plt.plot(np.arange(1,N+1),eG,"--",label="Gauss")
    plt.plot(np.arange(1,N+1),eP,label="Solver Python")
    plt.plot(np.arange(1,N+1),eL,label="LU")
    plt.xlabel("Nombre d'équations n")
    plt.ylabel("Erreur e_n")
    plt.legend()
    plt.title("Evolution de l'erreur en fonction du nombre d'équations")
    plt.show()
    return []
    
def TracerErreur2(N):
    eG=np.zeros(N)
    eP=np.zeros(N)
    eL=np.zeros(N)
    
    #Erreur
    for n in range(1,N+1):
        A=np.random.uniform(0,1,(n,n))#A est une matrice de nombres aléatoirement choisis entre 0 et 1
        while np.linalg.det(A)==0:
            A=np.random.uniform(0,1,(n,n))
        XT=np.ones((n,1))
        B=np.dot(A,XT)
        eG[n-1]=Norme(Gauss(A,B)-XT)
        eP[n-1]=Norme(np.linalg.solve(A,B)-XT)
        eL[n-1]=Norme(ResolLU(A,B)-XT)
    
    #Courbes
    plt.plot(np.arange(1,N+1),np.log(eG),"--",label="Gauss")
    plt.plot(np.arange(1,N+1),np.log(eP),label="Solver Python")
    plt.plot(np.arange(1,N+1),np.log(eL),label="LU")
    plt.xlabel("Nombre d'équations n")
    plt.ylabel("Erreur e_n")
    plt.legend()
    plt.title("Evolution de l'erreur en fonction du nombre d'équations")
    plt.show()
    return []

# Fonction pour la méthode de Jacobien (session5 Fonction  Complet)

def Jacobi(A,B,X0,n):
    #cette fonction calcule les n premières iterations, à partir de X0, de la suite de Jacobi pour le
    #système AX=B
    #Ces n+1 vecteurs Xi sont stockés en colonne dans la matrice Y 
    m=np.size(X0)
    Y=np.zeros((m,n+1))
    for i in range(m):
        Y[i,0]=X0[i]
    for k in range(n):
        for i in range(m):
            Y[i,k+1]=B[i]
            for j in range(m):
                if j != i:
                    Y[i,k+1]=Y[i,k+1]-A[i,j]*Y[j,k]
            Y[i,k+1]=Y[i,k+1]/A[i,i]
    return Y

####################################
    
def Norme(X):
	#X est un vecteur colonne ou un vecteur (3,)
	U=math.sqrt(np.dot(np.transpose(X),X)) #U est un vecteur (1,1)
	return U

####################################

def MatricesExemple(n):
    #cette fonction construit les matrices de l'exemple n-n du sujet
    A=np.zeros((n,n))
    for i in range(n-1):
        A[i,i]=3
        A[i,i+1]=-1
        A[i+1,i]=-2
    A[n-1,n-1]=3
    XT=np.ones((n,1))
    B=np.dot(A,XT)
    return (A,B,XT)
    
####################################

def TraceErreurJacobi(kmax):
    (A,B,XT)=MatricesExemple(10)
    X0=np.zeros((10,1))
    E=np.zeros(kmax)
    Y=Jacobi(A,B,X0,kmax)
    for k in range(kmax):
        E[k]=Norme(XT-Y[:,k].reshape(10,1))
    plt.plot(np.arange(kmax),np.log10(E),label="Méthode de Jacobi")
    plt.legend()
    plt.title="Evolution de l'erreur en fonction du nombre d'itérations"
    plt.show()
    
####################################

def GS(A,B,X0,n):
    m=np.size(X0)
    Y=np.zeros((m,n+1))
    for i in range(m):
        Y[i,0]=X0[i]
    k=0;
    while k<n:
        for i in range(m):
            Y[i,k+1]=B[i]
            for j in range(i):
                Y[i,k+1]=Y[i,k+1]-A[i,j]*Y[j,k+1]
            for j in range(i+1,m):
                Y[i,k+1]=Y[i,k+1]-A[i,j]*Y[j,k]
            Y[i,k+1]=Y[i,k+1]/A[i,i]
        k=k+1
    return Y

####################################    

def TraceErreurJ_GS(kmax):
    (A,B,XT)=MatricesExemple(10)
    X0=np.zeros((10,1))
    E=np.zeros(kmax)
    F=np.zeros(kmax)
    Y=Jacobi(A,B,X0,kmax)
    Z=GS(A,B,X0,kmax)
    for k in range(kmax):
        E[k]=Norme(XT-Y[:,k].reshape(10,1))
        F[k]=Norme(XT-Z[:,k].reshape(10,1))
    plt.plot(np.arange(kmax),np.log10(E),label="Méthode de Jacobi")
    plt.plot(np.arange(kmax),np.log10(F),label="Méthode de Gauss-Seidel")
    plt.legend()
    plt.title="Evolution de l'erreur en fonction du nombre d'itérations"
    plt.show()  