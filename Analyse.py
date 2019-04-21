from Corps import *
import time

   

TabTemp=[]
TabPres=[]

for loop in range(10):
    print(loop)
    #Nombre de corps
    nombrePlan=100
    
    #Durée de la simulation
    temps=0.1
    
    #Intervalle de temps (Il vaut mieux garder un multiple de 10 sinon, le ttab peut avoir des problemes de dimensionnement (N+1 colonnes plutot que N))
    dt=0.00001
    
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
    TailleInitiale=0.1
    VitesseInitiale=0+10*loop
    EcartType=0.5
    
    #Taille de la boite dans laquelle se passe les collisions (à garder STRICTEMENT inférieur à: TailleInitiale)
    TailleBoite=0.11

    #Generation des conditions initiales
    nombrePlan,TPosx,TPosy,TPosz,TVitx,TVity,TVitz=AttributionInitiale(TailleInitiale,VitesseInitiale,EcartType,nombrePlan,N,methode="Solide3D")
    
    
    ############################################################
    ###-----------------Programme Principale-----------------###
    ############################################################
    
    ProgrammePrincipal(TPosx,TPosy,TPosz,TVitx,TVity,TVitz,Epot,Ecin,Moment,TailleBoite,N,nombrePlan,dt)
    
    ttab2,Tpress,pression=Pression(Moment,TailleBoite,N,100,dt)
    
    temperature=Temperature(Ecin,N,nombrePlan)

    TabPres.append(pression)
    TabTemp.append(temperature)
TabPres=np.array(TabPres)
TabTemp=np.array(TabTemp)   
    
V=(2*TailleBoite)**3 
N=nombrePlan
coeff=np.polyfit(TabPres, TabTemp, 1)
print(coeff)
b=V/N-coeff[0]
a=coeff[1]/(N/V-b*(N/V)**2)

print('b= ',b)
print('a= ',a)    
TabTempTheo=TabPres*V/N    

plt.plot(TabPres,TabTemp,"o",label='Lennard-Jonnes calculé')
plt.plot(TabPres,TabTempTheo,label='Gaz parfait théorique')
plt.legend()
plt.show()



