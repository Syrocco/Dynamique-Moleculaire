from Corps import *
import time

#####################################################################
###-----------------Initialisation des paramètres-----------------###
#####################################################################


#Nombre de corps
nombrePlan=100

Tmasse=np.ones(nombrePlan)*0.0000000000001
Tmasse[0]=10**4

#Durée de la simulation
temps=10

#Intervalle de temps (Il vaut mieux garder un multiple de 10 sinon, le ttab peut avoir des problemes de dimensionnement (N+1 colonnes plutot que N))
dt=0.01

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
s=np.zeros((N,nombrePlan))
#Energie et quantité de mouvement
Epot=np.zeros(N)
Ecin=np.zeros(N)
Moment=np.zeros(N)  #Quantité de mouvement transmise aux parois, en valeur absolue puisqu'elle ne sert qu'à trouver une pression


#Demi longueur du cube dans lequel on place les corps et leurs vitesses + ecart type de la gaussienne en t=0
TailleInitiale=15
VitesseInitiale=0
EcartType=0


#Generation des conditions initiales
nombrePlan,TPosx,TPosy,TPosz,TVitx,TVity,TVitz=AttributionInitiale(TailleInitiale,VitesseInitiale,EcartType,nombrePlan,N,methode="Megalaxie",Masse=Tmasse[0])

"""
#Nombre de corps
nombrePlan=2500

Tmasse=np.ones(nombrePlan)*0.05
Tmasse[0]=10**4

#Durée de la simulation
temps=30

#Intervalle de temps (Il vaut mieux garder un multiple de 10 sinon, le ttab peut avoir des problemes de dimensionnement (N+1 colonnes plutot que N))
dt=0.003

temps=1898 secondes
"""
############################################################
###-----------------Programme Principale-----------------###
############################################################
t1=time.time()
ProgrammePrincipalGravite(TPosx,TPosy,TVitx,TVity,Epot,Ecin,Tmasse,N,nombrePlan,dt,s,Tmasse[1])
t2=time.time()
print(t2-t1)



#Calcul de l'énergie totale
Etot=Epot[1:-1]+Ecin[1:-1]







###########################################################
###-----------------Affichage Graphique-----------------###
###########################################################


if DispEne:
    plt.figure()
    plt.plot(ttab[:-1],Epot[1:-1],label="Epot")
    plt.plot(ttab[:-1],Ecin[1:-1],label="Ecin")
    plt.plot(ttab[:-1],Etot,label="Etot")
    plt.xlabel("Temps")
    plt.ylabel("Energie")
    plt.title("Graphe de l'evolution de l'énergie au cours du temps")
    plt.legend()
"""
fig = plt.figure()
axes = p3.Axes3D(fig)
"""


fig = plt.figure()
axes = plt.axes()
if Animation:
    ani = animation.FuncAnimation(fig, animate2D, fargs=(TPosx,TPosy,TailleInitiale,axes,s,Tmasse[1]), interval=100, save_count=int(N/10))
    if SaveAnimation:
        ani.save('./animationTest.mp4', fps=25,dpi=150)
else:
    for i in range(nombrePlan):
        axes.plot(TPosx[:,i],TPosy[:,i],TPosz[:,i])
        axes.plot(TPosx[:,0],TPosy[:,0],TPosz[:,0],"ro")


plt.show()
