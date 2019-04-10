
from Corps import *
import time

#####################################################################
###-----------------Initialisation des paramètres-----------------###
#####################################################################

#Nombre de corps
nombrePlan=10

#Durée de la simulation
temps=1

#Intervalle de temps (Il vaut mieux garder un multiple de 10 sinon, le ttab peut avoir des problemes de dimensionnement (N+1 colonnes plutot que N))
dt=0.005

#Nombre de simulation(s)
N=int(temps/dt)

#Options graphiques
DispEne=True
DispPression=True
DispMomentum=True
Animation=True
SaveAnimation=False



#Définition des tableaux contenant les différentes données
ttab=np.arange(dt,temps,dt)

#Energie et quantité de mouvement
Epot=np.zeros(N)
Ecin=np.zeros(N)
Moment=np.zeros(N)  #Quantité de mouvement transmise aux parois, en valeur absolue puisqu'elle ne sert qu'à trouver une pression


#Demi longueur du cube dans lequel on place les corps et leurs vitesses + ecart type de la gaussienne en t=0
TailleInitiale=3.5
VitesseInitiale=0
EcartType=4

#Taille de la boite dans laquelle se passe les collisions (à garder STRICTEMENT inférieur à: TailleInitiale)
TailleBoite=3

#Generation des conditions initiales
nombrePlan,TPosx,TPosy,TPosz,TVitx,TVity,TVitz=AttributionInitiale(TailleInitiale,VitesseInitiale,EcartType,nombrePlan,N,methode="Random")


############################################################
###-----------------Programme Principale-----------------###
############################################################

ProgrammePrincipal(TPosx,TPosy,TPosz,TVitx,TVity,TVitz,Epot,Ecin,Moment,TailleBoite,N,nombrePlan,dt)

ttab2,Tpress,pression=Pression(Moment,TailleBoite,N,100,dt)

temperature=Temperature(Ecin,N,nombrePlan)



#Calcul de l'énergie totale
Etot=Epot[1:-1]+Ecin[1:-1]

TVit = np.zeros((nombrePlan))
for i in range(nombrePlan):
    TVit[i] = norme(np.array([TVitx[-1,i], TVity[-1,i], TVitz[-1,i]]))
vitMoy=np.mean(TVit)
vit=np.linspace(0, 10 * vitMoy, 1000)

distrib=Maxwell(vit,temperature, vitMoy)
normalisateur=calculIntegral(Maxwell,0,10*vitMoy,10000,temperature, vitMoy)

plt.figure()
plt.hist(TVit, bins=int(nombrePlan/5),density=True)

plt.plot(vit,distrib/normalisateur)




###########################################################
###-----------------Affichage Graphique-----------------###
###########################################################
print("pression=", pression,"temperature=",temperature )

if DispEne:
    plt.figure()
    plt.plot(ttab[:-1],Epot[1:-1],label="Epot")
    plt.plot(ttab[:-1],Ecin[1:-1],label="Ecin")
    plt.plot(ttab[:-1],Etot,label="Etot")
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
    plt.plot(ttab[:-1],Moment[1:-1])
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
