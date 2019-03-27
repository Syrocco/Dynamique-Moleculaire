from Corps import *
import time

   

TabTemp=[]
TabPres=[]

for loop in range(50):
    print(loop)
    nombrePlan=30
    
    temps=150
    
    dt=0.01
    
    N=int(temps/dt)
    
    
    ttab=np.arange(dt,temps,dt)
    
    TPosx=np.zeros((N,nombrePlan))
    TPosy=np.zeros((N,nombrePlan))
    TPosz=np.zeros((N,nombrePlan))
    TVitx=np.zeros((N,nombrePlan))
    TVity=np.zeros((N,nombrePlan))
    TVitz=np.zeros((N,nombrePlan))
    
    Epot=np.zeros(N)
    Ecin=np.zeros(N)
    Moment=np.zeros(N) 
    TailleInitiale=2.5
    VitesseInitiale=0.01+0.1*loop
    EcartType=0
    
    TailleBoite=3
    
    AttributionInitiale(TailleInitiale,VitesseInitiale,EcartType,TPosx,TPosy,TPosz,TVitx,TVity,TVitz,nombrePlan)
    
      
              
    ProgrammePrincipal(TPosx,TPosy,TPosz,TVitx,TVity,TVitz,Epot,Ecin,Moment,TailleBoite,N,nombrePlan,dt)
    
    
    ttab2,Tpress,pression=Pression(Moment,TailleBoite,N,100,dt)
    temperature=Temperature(Ecin,N,nombrePlan)
    Etot=Epot[1:-1]+Ecin[1:-1]
    if abs(Etot[-2]-Etot[2])<5:
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



