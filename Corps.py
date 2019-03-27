import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rd
import mpl_toolkits.mplot3d.axes3d as p3
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import matplotlib.animation as animation
from numba import jit




################################################################
###-----------------Définition des fonctions-----------------###    
################################################################  
def cube(l,axes):
    points = np.array([[-l, -l, -l],
                      [l, -l, -l ],
                      [l, l, -l],
                      [-l, l, -l],
                      [-l, -l, l],
                      [l, -l, l],
                      [l, l, l],
                      [-l, l, l]])

    r = [-l,l]
    X, Y = np.meshgrid(r, r)    
    one =l* np.ones(4).reshape(2, 2)
    axes.plot_wireframe(X,Y,one, alpha=0.5)
    axes.plot_wireframe(X,Y,-one, alpha=0.5)
    axes.plot_wireframe(X,-one,Y, alpha=0.5)
    axes.plot_wireframe(X,one,Y, alpha=0.5)
    axes.plot_wireframe(one,X,Y, alpha=0.5)
    axes.plot_wireframe(-one,X,Y, alpha=0.5)
    axes.scatter3D(points[:, 0], points[:, 1], points[:, 2])    


def animate(i,PositionX,PositionY,PositionZ,demiLongueur,axes,l):
    i=i*10
    axes.clear()
    #axes.plot(PositionX[i,:1],PositionY[i,:1],PositionZ[i,:1],"ro")
    if i<300:
        axes.plot(PositionX[i,1:],PositionY[i,1:],PositionZ[i,1:],"ro")
    else:
        axes.plot(PositionX[i,1:],PositionY[i,1:],PositionZ[i,1:],"bo")    
    
    cube(l,axes)
    axes.set_xlim3d([-demiLongueur, demiLongueur])

    axes.set_ylim3d([-demiLongueur, demiLongueur])

    axes.set_zlim3d([-demiLongueur, demiLongueur])

def norme(Vec3D):
    return np.sqrt(normeCarree(Vec3D))

def normeCarree(Vec3D):
    return Vec3D[0]**2+Vec3D[1]**2+Vec3D[2]**2

#Avec la vitesse/température T et l'écart type E, la fonction renvoie UNE vitesse prise dans une gaussienne
def GenVitesse(T,E):
    vit_xyz=np.zeros(3)
    for i in range(3):
        vit_xyz[i]=rd.random()
    NormEtVit=norme(vit_xyz)/rd.normal(T,E)   #constante permettant de normer vit_xyz pour ensuite lui donner une norme prise au hasard dans une gaussienne

    for i in range(3):
        signe=rd.random()
        if signe<0.5:
            vit_xyz[i]=vit_xyz[i]/NormEtVit
        else:
            vit_xyz[i]=-vit_xyz[i]/NormEtVit
    
    return np.array(vit_xyz)

#La fonction renvoie la distribution des corps selon la demi-longueur envoyée
def GenPosition(nombre,rayon,methode):
    T=[]
    if methode=="Cube":
        Ecart=(2*rayon)/((nombre)**(1/3))
        for i in np.arange(-rayon,rayon,Ecart):
            for j in np.arange(-rayon+Ecart/2,rayon,Ecart):
                for z in np.arange(-rayon,rayon,Ecart):
                    T.append([i,j,z])
        if len(T)<nombre:
            print("L'attribution est mauvaise")
            return "L'attribution est mauvaise"
        rd.seed(1)
        return np.array(T[:nombre])


    if methode=="Random":
        for i in range(nombre):
            T.append(rd.random(3)*2*rayon-rayon)
    return np.array(T)


#La fonction attribue les positions et vitesses initiales aux corps                                                
def AttributionInitiale(rayon,Vitesse,Ecart,PositionX,PositionY,PositionZ,VitesseX,VitesseY,VitesseZ,nombre,methode="Cube"):
    particule=GenPosition(nombre, rayon, methode)
    for i in range(nombre):
        PositionX[0,i]=particule[i,0]
        PositionY[0,i]=particule[i,1]
        PositionZ[0,i]=particule[i,2]

    for i in range(nombre):
        A=GenVitesse(Vitesse,Ecart)
        VitesseX[0,i]=A[0]
        VitesseY[0,i]=A[1]
        VitesseZ[0,i]=A[2]
        


#Calcul de l'acceleration
@jit(nopython=True,cache=True)        
def CalculAcceleration(PositionX,PositionY,PositionZ,Corps,CorpsAutre,cpt):
	d=distance(PositionX[cpt-1,Corps],PositionX[cpt-1,CorpsAutre],PositionY[cpt-1,Corps],PositionY[cpt-1,CorpsAutre],PositionZ[cpt-1,Corps],PositionZ[cpt-1,CorpsAutre])
	a=((12/d**13)-(6/d**7))/1
	return a*(PositionX[cpt-1,Corps]-PositionX[cpt-1,CorpsAutre])/d, a*(PositionY[cpt-1,Corps]-PositionY[cpt-1,CorpsAutre])/d,a*(PositionZ[cpt-1,Corps]-PositionZ[cpt-1,CorpsAutre])/d, d 

#Calcul de la position et de la vitesse avec la méthode d'euler semi-explicite
@jit(nopython=True,cache=True)     
def CalculVitesseEtPosition(PositionX,PositionY,PositionZ,VitesseX,VitesseY,VitesseZ,AccelerationX,AccelerationY,AccelerationZ,Corps,cpt,dt):
    VitesseX[cpt,Corps]=VitesseX[cpt-1,Corps]+dt*AccelerationX
    VitesseY[cpt,Corps]=VitesseY[cpt-1,Corps]+dt*AccelerationY	
    VitesseZ[cpt,Corps]=VitesseZ[cpt-1,Corps]+dt*AccelerationZ
    PositionX[cpt,Corps]=PositionX[cpt-1,Corps]+dt*VitesseX[cpt,Corps]
    PositionY[cpt,Corps]=PositionY[cpt-1,Corps]+dt*VitesseY[cpt,Corps]
    PositionZ[cpt,Corps]=PositionZ[cpt-1,Corps]+dt*VitesseZ[cpt,Corps] 
    
    
#Teste si la particule à l'étape i se trouve dans la boite ou non                
def DansBoite(demiLongueur,PositionX,PositionY,PositionZ, Corps, cpt):
        
    #r la demi longueur du cube
    if PositionX[cpt,Corps]<-demiLongueur:
        return "x-"
    if PositionX[cpt,Corps]>demiLongueur:
        return "x+"
    if PositionY[cpt,Corps]<-demiLongueur:
        return "y-"
    if PositionY[cpt,Corps]>demiLongueur:
        return "y+"
    if PositionZ[cpt,Corps]<-demiLongueur:
        return "z-"
    if PositionZ[cpt,Corps]>demiLongueur:
        return "z+"
    return("no")

#Modifie, selon le resultat de la fonction DansBoite(), la vitesse et la position de la particule à l'étape i pour simuler une collision, enregistre aussi la quantité de mouvement transmise aux parois
@jit(nopython=True,cache=True) 
def modif(info,demiLongueur,PositionX,PositionY,PositionZ,VitesseX,VitesseY,VitesseZ,QuantDeMouv,Corps,cpt):
    if info=="no":
        return()
    if info=='x-':
        PositionX[cpt,Corps]=PositionX[cpt,Corps] + 2*(-demiLongueur-PositionX[cpt,Corps])
        QuantDeMouv[cpt]=QuantDeMouv[cpt]+abs(2*VitesseX[cpt,Corps])
        VitesseX[cpt,Corps]=-VitesseX[cpt,Corps]
        return()
    if info=='x+':
        PositionX[cpt,Corps]=PositionX[cpt,Corps] + 2*(demiLongueur-PositionX[cpt,Corps])
        QuantDeMouv[cpt]=QuantDeMouv[cpt]+abs(2*VitesseX[cpt,Corps])
        VitesseX[cpt,Corps]=-VitesseX[cpt,Corps]
        return()
    if info=='y-':
        PositionY[cpt,Corps]=PositionY[cpt,Corps] + 2*(-demiLongueur-PositionY[cpt,Corps])
        QuantDeMouv[cpt]=QuantDeMouv[cpt]+abs(2*VitesseY[cpt,Corps])
        VitesseY[cpt,Corps]=-VitesseY[cpt,Corps]
        return()
    if info=='y+':
        PositionY[cpt,Corps]=PositionY[cpt,Corps] + 2*(demiLongueur-PositionY[cpt,Corps])
        QuantDeMouv[cpt]=QuantDeMouv[cpt]+abs(2*VitesseY[cpt,Corps])
        VitesseY[cpt,Corps]=-VitesseY[cpt,Corps]
        return()
    if info=='z-':
        PositionZ[cpt,Corps]=PositionZ[cpt,Corps] + 2*(-demiLongueur-PositionZ[cpt,Corps])
        QuantDeMouv[cpt]=QuantDeMouv[cpt]+abs(2*VitesseZ[cpt,Corps])
        VitesseZ[cpt,Corps]=-VitesseZ[cpt,Corps]
        return()
    if info=='z+':
        PositionZ[cpt,Corps]=PositionZ[cpt,Corps] + 2*(demiLongueur-PositionZ[cpt,Corps])
        QuantDeMouv[cpt]=QuantDeMouv[cpt]+abs(2*VitesseZ[cpt,Corps])
        VitesseZ[cpt,Corps]=-VitesseZ[cpt,Corps]
        
@jit(nopython=True,cache=True) 
def distance(PositionX1,PositionX2,PositionY1,PositionY2,PositionZ1,PositionZ2):
	return np.sqrt((PositionX1-PositionX2)**2+(PositionY1-PositionY2)**2+(PositionZ1-PositionZ2)**2)

@jit(nopython=True,cache=True) 
def Ecinetique(VitesseX,VitesseY,VitesseZ):
	return 0.5*1*(VitesseX**2+VitesseY**2+VitesseZ**2)

def Epotentielle(distance):
	return 0.5*((1/distance**12)-(1/distance**6))   #Multiplication par 0.5 parceque je compte deux fois l'energie: pour un couple {i,j} de particule, je compte Eij et Eji, donc il faut diviser par deux

#Calcul de la pression exercée sur les parois 
def Pression(QuantDeMouv,demiLongueur,nombreDiteration,pas,dt):
    aireBoite=6*(2*demiLongueur)**2
    nombreDivision=int(nombreDiteration/pas)

    ttab2=np.linspace(0,dt*nombreDivision,nombreDivision)
    TMoment=np.zeros(nombreDivision)
    for i in range(nombreDivision):
        TMoment[i]=sum(QuantDeMouv[i*pas:(i*pas)+pas])
    pression=TMoment/(pas*dt*aireBoite)
    pression=np.array(pression)
    pressionMoyenne=np.mean(pression[int(len(pression)*9/10):])
    return ttab2,pression,pressionMoyenne

#Calcul temperature
def Temperature(EnergieCinetique,nombreDiteration,nombre):
    kb=1
    eneCinMoyenne=0
    for i in range(int(nombreDiteration/2),nombreDiteration):
        eneCinMoyenne+=2*EnergieCinetique[i]/nombreDiteration
    return eneCinMoyenne/(1.5*kb*nombre)


def ProgrammePrincipal(PositionX,PositionY,PositionZ,VitesseX,VitesseY,VitesseZ,EnergiePotentielle,EnergieCinetique,QuantDeMouv,demiLongueur,nombreDiteration,nombre,dt,EnergieVoulue=0,coeff=1.001):
    print("Début des calculs")
    for i in range(1,nombreDiteration):
        #print(i,"/",nombreDiteration)
        for Corps in range(nombre):
            ax=0
            ay=0
            az=0
            for CorpsAutre in range(nombre):
               if CorpsAutre!=Corps:      
                   a=CalculAcceleration(PositionX,PositionY,PositionZ,Corps,CorpsAutre,i)
                   ax+=a[0]
                   ay+=a[1]
                   az+=a[2]
                   EnergiePotentielle[i-1]=EnergiePotentielle[i-1]+Epotentielle(a[3])
            CalculVitesseEtPosition(PositionX,PositionY,PositionZ,VitesseX,VitesseY,VitesseZ,ax,ay,az,Corps,i,dt)
            modif(DansBoite(demiLongueur,PositionX,PositionY,PositionZ,Corps, i),demiLongueur,PositionX,PositionY,PositionZ,VitesseX,VitesseY,VitesseZ,QuantDeMouv,Corps,i)    #A mettre en commentaire pour desactiver les collisions
            EnergieCinetique[i]=EnergieCinetique[i]+Ecinetique(VitesseX[i,Corps],VitesseY[i,Corps],VitesseZ[i,Corps])

            N=10
            if i%N==0 and i>=N and 300<i<10000:

                EnergieMoyenne=np.mean(EnergieCinetique[i-N:i])
                if EnergieVoulue-EnergieMoyenne>0.01:
                    for n in range(nombre):
                        VitesseX[i,n]=VitesseX[i,n]*coeff
                        VitesseY[i,n]=VitesseY[i,n]*coeff
                        VitesseZ[i,n]=VitesseZ[i,n]*coeff
                elif EnergieVoulue-EnergieMoyenne<0.01:
                    for n in range(nombre):
                        VitesseX[i,n]=VitesseX[i,n]/coeff
                        VitesseY[i,n]=VitesseY[i,n]/coeff
                        VitesseZ[i,n]=VitesseZ[i,n]/coeff 
