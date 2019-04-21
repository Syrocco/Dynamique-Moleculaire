
from Corps import *
import time

#############################################################################################
###-----------------Initialisation des paramètres et Programme Principal-----------------###
#############################################################################################

#Nombre de corps
nombrePlan=100

#Durée de la simulation
temps=10

dt=0.01
#Nombre de simulation(s)
N=int(temps/dt)

#Options graphiques et sauvegarde
DispEne=True
DispPression=True
DispMomentum=True
DispDistributionVit=True
Animation=True
SaveAnimation=False
SaveTab=True



#Définition des tableaux contenant les différentes données
ttab=np.linspace(dt,temps,N)

#Energie et quantité de mouvement
Epot=np.zeros(N)
Ecin=np.zeros(N)
Moment=np.zeros(N)  #Quantité de mouvement transmise aux parois, en valeur absolue puisqu'elle ne sert qu'à trouver une pression


#Demi longueur du cube dans lequel on place les corps et leurs vitesses + ecart type de la gaussienne en t=0
TailleInitiale=4
VitesseInitiale=5
EcartType=0

#Taille de la boite dans laquelle se passe les collisions (à garder STRICTEMENT inférieur à: TailleInitiale)
TailleBoite=5

#Generation des conditions initiales
nombrePlan,TPosx,TPosy,TPosz,TVitx,TVity,TVitz=AttributionInitiale(TailleInitiale,VitesseInitiale,EcartType,nombrePlan,N,methode="Cube")



ProgrammePrincipal(TPosx,TPosy,TPosz,TVitx,TVity,TVitz,Epot,Ecin,Moment,TailleBoite,N,nombrePlan,dt)


ttab2,Tpress,pression=Pression(Moment,TailleBoite,N,100,dt)
temperature=Temperature(Ecin,N,nombrePlan)
Etot=Epot[1:-1]+Ecin[1:-1]


#########################################################################
###-----------------Affichage Graphique et sauvegarde-----------------###
#########################################################################
print("pression=", pression,"temperature=",temperature )

if SaveTab:
    Donnees=np.zeros(N-2,dtype=[('ttab',float),('Epot',float),('Ecin',float),('Moment',float)])
    Donnees['ttab']=ttab[1:-1]
    Donnees['Epot']=Epot[1:-1]
    Donnees['Ecin']=Ecin[1:-1]
    Donnees['Moment']=Moment[1:-1]
    information='Taille initiale='+str(TailleInitiale)+'    Vitesse initiale='+str(VitesseInitiale)+' et '+str(EcartType)+'    Taille boîte='+str(TailleBoite)+'    dt='+str(dt)+'    nombre='+str(nombrePlan)
    nom='nb='+str(nombrePlan)+',vit='+str(VitesseInitiale)+',TailleBoite='+str(TailleBoite)+'.txt'
    np.savetxt(nom,Donnees,header=information)
   
    


if DispDistributionVit:
    TVit = np.zeros((nombrePlan))
    for i in range(nombrePlan):
        TVit[i] = norme(np.array([TVitx[-1,i], TVity[-1,i], TVitz[-1,i]]))
    vitMoy=np.mean(TVit)
    vit=np.linspace(0, 3 * vitMoy, 1000)
    
    distrib=Maxwell(vit,temperature)

    plt.figure()
    plt.hist(TVit, bins=int(nombrePlan/10),density=True)
    
    plt.plot(vit,distrib)


if DispEne:
    plt.figure()
    plt.plot(ttab[1:-1],Epot[1:-1],label="Epot")
    plt.plot(ttab[1:-1],Ecin[1:-1],label="Ecin")
    plt.plot(ttab[1:-1],Etot,label="Etot")
    plt.xlabel("Temps")
    plt.ylabel("Energie")
    plt.title("Graphe de l'evolution de l'énergie au cours du temps")
    plt.legend()

if DispPression:
    plt.figure()
    plt.plot(ttab2,Tpress)
    plt.xlabel("Temps")
    plt.ylabel("Pression")
    plt.title("Graphe de l'evolution de la pression dans la boîte au cours du temps")
if DispMomentum:
    plt.figure()
    plt.plot(ttab[1:-1],Moment[1:-1])
    plt.xlabel("Temps")
    plt.ylabel("Quantité de mouvement")
    plt.title("Graphe de l'evolution de la pression dans la boîte au cours du temps")


if Animation:
    fig = plt.figure()
    axes = p3.Axes3D(fig)
    ani = animation.FuncAnimation(fig, animate, fargs=(TPosx,TPosy,TPosz,TailleBoite,axes), interval=100, save_count=int(N/10))
    if SaveAnimation:
        ani.save('./animationTrash.mp4', fps=40, dpi=150)
else:
    fig = plt.figure()
    axes = p3.Axes3D(fig)
    cube(TailleBoite,axes)
    for i in range(nombrePlan):
        axes.plot(TPosx[:,i],TPosy[:,i],TPosz[:,i],"o")





plt.show()
