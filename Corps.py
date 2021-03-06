import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rd
import mpl_toolkits.mplot3d.axes3d as p3
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import matplotlib.animation as animation
from numba import jit





###########################################################
###-----------------Fonctions: Fourre tout--------------###
###########################################################
def norme(Vec3D):
    return np.sqrt(normeCarree(Vec3D))

def norme2D(Vec2D):
    return np.sqrt(normeCarree2D(Vec2D))


def normeCarree(Vec3D):
    return Vec3D[0]**2+Vec3D[1]**2+Vec3D[2]**2

def normeCarree2D(Vec2D):
    return Vec2D[0]**2+Vec2D[1]**2

@jit(nopython=True,cache=True)
def distance(PositionX1,PositionX2,PositionY1,PositionY2,PositionZ1,PositionZ2):
	return np.sqrt((PositionX1-PositionX2)**2+(PositionY1-PositionY2)**2+(PositionZ1-PositionZ2)**2)


@jit(nopython=True,cache=True)
def distance2D(PositionX1,PositionX2,PositionY1,PositionY2):
	return np.sqrt((PositionX1-PositionX2)**2+(PositionY1-PositionY2)**2)



####################################################################
###-----------------Fonctions: Conditions initiales--------------###
####################################################################



#Avec la vitesse mu et l'écart type sigma, la fonction renvoie UNE vitesse prise dans une gaussienne
def GenVitesse(mu,sigma):
    vit_xyz=np.zeros(3)
    if (mu == 0 and sigma == 0):
        return np.array(vit_xyz)

    for i in range(3):
        vit_xyz[i]=rd.random()
    NormEtVit=norme(vit_xyz)/rd.normal(mu,sigma)   #constante permettant de normer vit_xyz pour ensuite lui donner une norme prise au hasard dans une gaussienne

    for i in range(3):
        signe=rd.random()
        if signe<0.5:
            vit_xyz[i]=vit_xyz[i]/NormEtVit
        else:
            vit_xyz[i]=-vit_xyz[i]/NormEtVit

    return np.array(vit_xyz)





#La fonction renvoie la distribution des corps 
def GenPosition(nombre,rayon,methode):
    T=[]

    if methode=="Cube":
        Ecart=(2*rayon)/((nombre)**(1/3))
        for i in np.arange(-rayon,rayon,Ecart):
            for j in np.arange(-rayon+Ecart/2,rayon,Ecart):
                for z in np.arange(-rayon,rayon,Ecart):
                    T.append([i,j,z])
        if len(T)<nombre:
            return "L'attribution est mauvaise"
        #rd.seed(1)
        return np.array(T[:nombre])

   
    if methode=="Random":
        T.append(rd.random(3)*2*rayon-rayon)
        while len(T)!=nombre:
            a=rd.random(3)*2*rayon-rayon
            for i in range(len(T)):             
                if distance(np.array(T)[i,0],a[0],np.array(T)[i,1],a[1],np.array(T)[i,2],a[2])<(2**(1/6))/1.3:
                    break
                if i==len(T)-1:
                    T.append(a)
        return np.array(T)
            


    if methode=="Solide2D":
        dist=2**(1/6)
        alt=0
        for i in np.arange(-rayon,rayon,dist):
                for j in np.arange(-rayon+alt*dist/2,rayon+alt*dist/2,dist):
                        T.append([i,j,0])
                if alt==1:
                    alt=0
                else:
                    alt=1

        return np.array(T), len(T)

    if methode=="Solide3D":

        dist=2**(1/6)

        for z in np.arange(-rayon,rayon,2*dist):
            alt=0
            for i in np.arange(-rayon,rayon,dist):
                for j in np.arange(-rayon+alt*dist/2,rayon+alt*dist/2,dist):
                        T.append([i,j,z])
                if alt==1:
                    alt=0
                else:
                    alt=1
            alt=1
            for i in np.arange(-rayon,rayon,dist):
                for j in np.arange(-rayon+alt*dist/2,rayon+alt*dist/2,dist):
                        T.append([i,j,z+dist])
                if alt==1:
                    alt=0
                else:
                    alt=1

        return np.array(T), len(T)


    if methode=="Megalaxie":
        T.append([0,0,0])
        while(len(T)!=nombre):
            x=rd.random(2) * (2*rayon)-rayon
            if (rayon/7 <= norme2D(x) <= rayon):
                T.append([x[0],x[1],0])
        return np.array(T)


#La fonction attribue les positions et vitesses initiales aux corps
def AttributionInitiale(rayon,Vitesse,Ecart,nombre,N,methode="Cube",G=1,Masse=1):
    if methode=="Solide3D" or methode=="Solide2D":
        particule, nombre=GenPosition(nombre, rayon, methode)
    else:
        particule=GenPosition(nombre, rayon, methode)

    #Position et vitesse
    PositionX=np.zeros((N,nombre))
    PositionY=np.zeros((N,nombre))
    PositionZ=np.zeros((N,nombre))
    VitesseX=np.zeros((N,nombre))
    VitesseY=np.zeros((N,nombre))
    VitesseZ=np.zeros((N,nombre))

    for i in range(nombre):
        PositionX[0,i]=particule[i,0]
        PositionY[0,i]=particule[i,1]
        PositionZ[0,i]=particule[i,2]

    if methode!="Megalaxie":
        for i in range(1,nombre):
            A=GenVitesse(Vitesse,Ecart)
            VitesseX[0,i]=A[0]
            VitesseY[0,i]=A[1]
            VitesseZ[0,i]=A[2]

    else:
        for i in range(1,nombre):
            vector=np.array([PositionX[0,i],PositionY[0,i],0])
            Vnormal=np.cross(vector, np.array([0, 0, 1]))
            Vnormal=Vnormal/norme(Vnormal)

            #print(vector,Vnormal)
            VTot=np.sqrt((G*Masse)/norme(vector))
            VitesseX[0,i]=Vnormal[0]*VTot
            VitesseY[0,i]=Vnormal[1]*VTot
            VitesseZ[0,i]=Vnormal[2]*VTot


    return nombre,PositionX,PositionY,PositionZ,VitesseX,VitesseY,VitesseZ



#############################################################################
###-----------------Fonctions: Calculs principaux et mesures--------------###
#############################################################################

#-------------Resolution d'equa diff-------------#

#Calcul de l'acceleration
@jit(nopython=True,cache=True)
def CalculAcceleration(PositionX,PositionY,PositionZ,Corps,CorpsAutre,cpt):
    d=distance(PositionX[cpt-1,Corps],PositionX[cpt-1,CorpsAutre],PositionY[cpt-1,Corps],PositionY[cpt-1,CorpsAutre],PositionZ[cpt-1,Corps],PositionZ[cpt-1,CorpsAutre])
    a=((12/d**13)-(6/d**7))
    return a*(PositionX[cpt-1,Corps]-PositionX[cpt-1,CorpsAutre])/d, a*(PositionY[cpt-1,Corps]-PositionY[cpt-1,CorpsAutre])/d,a*(PositionZ[cpt-1,Corps]-PositionZ[cpt-1,CorpsAutre])/d, d

@jit(nopython=True,cache=True)
def CalculAccelerationGravite(PositionX,PositionY,Tmasse,Corps,CorpsAutre,cpt):
    d=distance2D(PositionX[cpt-1,Corps],PositionX[cpt-1,CorpsAutre],PositionY[cpt-1,Corps],PositionY[cpt-1,CorpsAutre])
    a=-Tmasse[CorpsAutre]/d**2
    return a*(PositionX[cpt-1,Corps]-PositionX[cpt-1,CorpsAutre])/d, a*(PositionY[cpt-1,Corps]-PositionY[cpt-1,CorpsAutre])/d, d




#Calcul de la position et de la vitesse avec la méthode d'euler semi-explicite
@jit(nopython=True,cache=True)
def CalculVitesseEtPosition(PositionX,PositionY,PositionZ,VitesseX,VitesseY,VitesseZ,AccelerationX,AccelerationY,AccelerationZ,Corps,cpt,dt):
    VitesseX[cpt,Corps]=VitesseX[cpt-1,Corps]+dt*AccelerationX
    VitesseY[cpt,Corps]=VitesseY[cpt-1,Corps]+dt*AccelerationY
    VitesseZ[cpt,Corps]=VitesseZ[cpt-1,Corps]+dt*AccelerationZ
    PositionX[cpt,Corps]=PositionX[cpt-1,Corps]+dt*VitesseX[cpt,Corps]
    PositionY[cpt,Corps]=PositionY[cpt-1,Corps]+dt*VitesseY[cpt,Corps]
    PositionZ[cpt,Corps]=PositionZ[cpt-1,Corps]+dt*VitesseZ[cpt,Corps]

@jit(nopython=True,cache=True)
def CalculVitesseEtPositionGravite(PositionX,PositionY,VitesseX,VitesseY,AccelerationX,AccelerationY,Corps,cpt,dt):
    VitesseX[cpt,Corps]=VitesseX[cpt-1,Corps]+dt*AccelerationX
    VitesseY[cpt,Corps]=VitesseY[cpt-1,Corps]+dt*AccelerationY
    PositionX[cpt,Corps]=PositionX[cpt-1,Corps]+dt*VitesseX[cpt,Corps]
    PositionY[cpt,Corps]=PositionY[cpt-1,Corps]+dt*VitesseY[cpt,Corps]


#Calcul de la position et de la vitesse au temps t1 pour utiliser ensuite verlet
def PosVit_t1(PositionX,PositionY,PositionZ,VitesseX,VitesseY,VitesseZ,nombre,dt):
    for planete in range(nombre):
        ax=0
        ay=0
        az=0
        for planeteAutre in range(nombre):
            if planeteAutre!=planete:
                d=distance(PositionX[0,planete],PositionX[0,planeteAutre],PositionY[0,planete],PositionY[0,planeteAutre],PositionZ[0,planete],PositionZ[0,planeteAutre])
                a=((12/d**13)-(6/d**7))
                ax+=a*(PositionX[0,planete]-PositionX[0,planeteAutre])/d
                ay+=a*(PositionY[0,planete]-PositionY[0,planeteAutre])/d
                az+=a*(PositionZ[0,planete]-PositionZ[0,planeteAutre])/d


        VitesseX[1,planete]=VitesseX[0,planete]+dt*ax
        VitesseY[1,planete]=VitesseY[0,planete]+dt*ay
        VitesseZ[1,planete]=VitesseZ[0,planete]+dt*ax
        PositionX[1,planete]=PositionX[0,planete]+dt*VitesseX[1,planete]
        PositionY[1,planete]=PositionY[0,planete]+dt*VitesseY[1,planete]
        PositionZ[1,planete]=PositionZ[0,planete]+dt*VitesseZ[1,planete]

#Calcul de la position et de la vitesse avec la méthode de verlet (pas spécialement adaptée pour les collisions, je ne l'utilise pas)
@jit(nopython=True,cache=True)
def VerletPosition(PositionX,PositionY,PositionZ,VitesseX,VitesseY,VitesseZ,AccelerationX,AccelerationY,AccelerationZ,Corps,cpt,dt):
    PositionX[cpt+1,Corps]=2*PositionX[cpt,Corps]-PositionX[cpt-1,Corps]+AccelerationX*dt**2
    PositionY[cpt+1,Corps]=2*PositionY[cpt,Corps]-PositionY[cpt-1,Corps]+AccelerationY*dt**2
    PositionZ[cpt+1,Corps]=2*PositionZ[cpt,Corps]-PositionZ[cpt-1,Corps]+AccelerationZ*dt**2
    VitesseX[cpt,Corps]=(PositionX[cpt+1,Corps]-PositionX[cpt-1,Corps])/(2*dt)
    VitesseY[cpt,Corps]=(PositionY[cpt+1,Corps]-PositionY[cpt-1,Corps])/(2*dt)
    VitesseZ[cpt,Corps]=(PositionZ[cpt+1,Corps]-PositionZ[cpt-1,Corps])/(2*dt)




#-------------Tests divers, mesures et modification des données-------------#


#Teste si la particule à l'étape i se trouve dans la boite ou non
@jit(nopython=True,cache=True)
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
def Ecinetique(VitesseX,VitesseY,VitesseZ):
	return 0.5*1*(VitesseX**2+VitesseY**2+VitesseZ**2)

@jit(nopython=True,cache=True)
def EcinetiqueGravite(VitesseX,VitesseY,m):
	return 0.5*m*(VitesseX**2+VitesseY**2)

@jit(nopython=True,cache=True)
def Epotentielle(distance):
	return 0.5*((1/distance**12)-(1/distance**6))

@jit(nopython=True,cache=True)
def EpotentielleGravite(distance,m1,m2):
    return -0.5*m1*m2/distance


@jit(nopython=True,cache=True)
def TempModif(EnergieCinetique,VitesseX,VitesseY,VitesseZ,i,nombre,nombreDiteration):
    N=10
    EnergieVoulue=1
    coeff=1.001
    if i%N==0 and i>=2*N and i<9*nombreDiteration/10:
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

@jit(nopython=True,cache=True)
def collisions(PositionX,PositionY,VitesseX,VitesseY,Tmasse,cpt,MasseNormee):   #I am become death, the destroyer of optimization :^)
    for Corps in range(len(Tmasse)):
        for CorpsAutre in range(len(Tmasse)):
            if Corps!=CorpsAutre:
                if distance2D(PositionX[cpt,Corps],PositionX[cpt,CorpsAutre],PositionY[cpt,Corps],PositionY[cpt,CorpsAutre])<=0.5:
                    VitesseX[cpt,Corps]=(Tmasse[Corps]*VitesseX[cpt,Corps]+Tmasse[CorpsAutre]*VitesseX[cpt,CorpsAutre])/(Tmasse[Corps]+Tmasse[CorpsAutre])
                    VitesseY[cpt,Corps]=(Tmasse[Corps]*VitesseY[cpt,Corps]+Tmasse[CorpsAutre]*VitesseY[cpt,CorpsAutre])/(Tmasse[Corps]+Tmasse[CorpsAutre])
                    Tmasse[Corps]=Tmasse[Corps]+Tmasse[CorpsAutre]
                    PositionX[cpt,CorpsAutre]=500*(1+cpt)+100*(1+rd.random())
                    PositionY[cpt,CorpsAutre]=500*(1+cpt)+100*(1+rd.random())
                    VitesseX[cpt,CorpsAutre]=0
                    VitesseY[cpt,CorpsAutre]=0


@jit(nopython=True,cache=True)
def ProgrammePrincipal(PositionX,PositionY,PositionZ,VitesseX,VitesseY,VitesseZ,EnergiePotentielle,EnergieCinetique,QuantDeMouv,demiLongueur,nombreDiteration,nombre,dt):
    print("Début des calculs")
    for i in range(1,nombreDiteration):
        print(i,"/",nombreDiteration)
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
            #TempModif(EnergieCinetique,VitesseX,VitesseY,VitesseZ,i,nombre,nombreDiteration)


@jit(nopython=True,cache=True)
def ProgrammePrincipalGravite(PositionX,PositionY,VitesseX,VitesseY,EnergiePotentielle,EnergieCinetique,Tmasse,nombreDiteration,nombre,dt,s,MasseNormee):
    print("Début des calculs")
    for i in range(1,nombreDiteration):
        print(i,"/",nombreDiteration)
        Corps=0
        for Corps in range(nombre):
            ax=0
            ay=0
            CorpsAutre=0
            for CorpsAutre in range(nombre):
                if CorpsAutre!=Corps:
                    a=CalculAccelerationGravite(PositionX,PositionY,Tmasse,Corps,CorpsAutre,i)
                    ax+=a[0]
                    ay+=a[1]
                    EnergiePotentielle[i-1]=EnergiePotentielle[i-1]+EpotentielleGravite(a[2],Tmasse[Corps],Tmasse[CorpsAutre])
            CalculVitesseEtPositionGravite(PositionX,PositionY,VitesseX,VitesseY,ax,ay,Corps,i,dt)
            EnergieCinetique[i]=EnergieCinetique[i]+EcinetiqueGravite(VitesseX[i,Corps],VitesseY[i,Corps],Tmasse[Corps])
        collisions(PositionX,PositionY,VitesseX,VitesseY,Tmasse,i,MasseNormee)
        s[i,:]=Tmasse


############################################################
###-----------------Fonctions: Post calculs--------------###
############################################################

#-------------Utilisation des données-------------#

#Calcul temperature
def Temperature(EnergieCinetique,nombreDiteration,nombre):
    kb=1
    eneCinMoyenne=0
    for i in range(int(nombreDiteration/2),nombreDiteration):
        eneCinMoyenne+=2*EnergieCinetique[i]/nombreDiteration
    return eneCinMoyenne/(1.5*kb*nombre)


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



def Maxwell(x,temp):
    return np.sqrt(2/np.pi)*(1/temp)**(3/2)*x**2* np.exp(-(x)**2/(2*temp))


#-------------Affichage graphique-------------#


def animate(i,PositionX,PositionY,PositionZ,demiLongueur,axes):
    i=i*10
    axes.clear()
    axes.view_init(15, i/30)
    axes.scatter(PositionX[i,1:],PositionY[i,1:],PositionZ[i,1:],"bo")
    #axes.scatter(PositionX[i,:1],PositionY[i,:1],PositionZ[i,:1],"ro")


    cube(demiLongueur,axes)
    axes.set_xlim3d([-demiLongueur, demiLongueur])

    axes.set_ylim3d([-demiLongueur, demiLongueur])

    axes.set_zlim3d([-demiLongueur, demiLongueur])

def animate2D(i,PositionX,PositionY,demiLongueur,axes,s,MasseNormee):
    i=i*10
    axes.clear()
    axes.scatter(PositionX[i,:1],PositionY[i,:1],s=250,c="#FF0000")
    axes.scatter(PositionX[i,1:],PositionY[i,1:],s=0.5/MasseNormee*s[i,1:],alpha=0.8,c="#4B0082")


    axes.axis([-demiLongueur, demiLongueur, -demiLongueur, demiLongueur])



def cube(l,axes):     #fonctions prise ici: https://stackoverflow.com/questions/33540109/plot-surfaces-on-a-cube/33542678
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
    axes.plot_wireframe(X,Y,one)
    axes.plot_wireframe(X,Y,-one)
    axes.plot_wireframe(X,-one,Y)
    axes.plot_wireframe(X,one,Y)
    axes.plot_wireframe(one,X,Y)
    axes.plot_wireframe(-one,X,Y)
    axes.scatter3D(points[:, 0], points[:, 1], points[:, 2])
