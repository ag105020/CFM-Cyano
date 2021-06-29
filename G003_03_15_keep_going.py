'''
Created on Apr 8, 2019


Here things are normalized by QcBio: biomass C
@author: keiin
'''

from pylab import *
from FigSetting2 import *
from Savefig3 import *
from Savetxt import *
rcParams['axes.prop_cycle'] = cycler('color',['b', 'r', 'g'])
rcParams.update({'font.size': 15})
rcParams.update({'lines.linewidth': 2.5})
rcParams.update({'figure.autolayout': True})
rcParams['figure.figsize']=5,4
rcParams.update({'axes.linewidth':1.2}) 
rcParams.update({'xtick.major.width':1.2})  #
rcParams.update({'ytick.major.width':1.2})  #
rcParams.update({'lines.markeredgewidth': 1}) #default = 1
rcParams.update({'ytick.right': 'True'})
rcParams.update({'ytick.direction': 'in'})
rcParams.update({'xtick.top': 'True'})
rcParams.update({'xtick.direction': 'in'})
rcParams['xtick.major.pad']='8'

#=====================================
# Defining parameters
#=====================================

CO2air = 1e6   #(uM) DIC concentration when equilibrium with air

Eno3 = 0.774  #From a800 05 12 07 Phytoplankton model 
Enstore = 0.383 #From 317 12 05 with epNH4=0.6 and C5H7O2N1
E = Eno3

#----------------------
# Storage
#----------------------
pCstoreMax = 8 #(d-1)
CstoreMax = 0.65  #(molC molC-1)

#:::::::::::::::::::::::::::::::::::::::::::::::
# Key parameters; These are really keys
#:::::::::::::::::::::::::::::::::::::::::::::::
PhoMax = 4.5    #(d-1) Maximum photosynthetsis rate
Kco2 = 250     #(uM)

Hco2 = 1.65 #(d-1) gas exchange constant 
Cvco2 = 2.9e3 #(dimensionless) factor for adjusting CO2 consumption

Cstoreloss = 0.6 #(d-1) Cstoreloss due to its conversion to N containing molecules (biomass)

#============================
# 
#============================
dts = 10  #(s) dt in second
dt = dts/86400 #(d) dt in day
t = arange(0,86400/24*24+dts,dts)/86400 #(d) time
th = t*24 #(h) time
U = arange(size(t))
NightTime = 3600*14/86400
#============================

QcBio = 1

#===============================
# Preparing for arrays 
#===============================
def upd1(f,df):  #upd = update
    f[i+1] = f[i] + df*dt 
    return f

def o(): 
    return zeros((size(t)))

Pho = o()
CO2 = o()
Cstore = o()
Cstore2 = o()
Qc = o()
X = o() #(uM C) proxy of population (assuming constant cell C of QcBio).
Mu = o()
MuC = o()
pCstore = o()

#=====================
# Initial conditions
#=====================
CO2[0] = 812.06 #(uM)
CO2 = CO2air*ones(size(o()))
Cstore[0] = 0.38 #(molC molC-1)
Cstore2[0] = Cstore[0]
X[0] = 1

#=====================

for i in U[:-1]:
    if t[i] < NightTime:
        Pho[i] = PhoMax*CO2[i]/(CO2[i]+Kco2)
        pCstore[i] = pCstoreMax*(CstoreMax - Cstore[i])/CstoreMax
        Mu[i] = (Pho[i] - pCstore[i]*(1+E))/(1+(1-Cstore[i])*E)
    
    if Pho[i] <= 0:
        Pho[i] = 0
        pCstore[i] = 0
        Mu[i] = 0
    
    if Pho[i] < pCstore[i]*(1+E):
        pCstore[i] = Pho[i]/(1+E)
        Mu[i] = 0
    
    Qc[i] = QcBio + Cstore[i]
    
    if t[i]<NightTime:
       # dCstore = Pho[i] - Mu[i]*(1+E)-Mu[i]*(Qc[i]-1)
        dCstore = pCstore[i] - Mu[i]*Cstore[i]
        dCstore2 = dCstore
      #  dCstore = pCstore[i]
    
        dX = Mu[i]*X[i]#/Qc[i]
        dCO2 = Hco2*(CO2air - CO2[i]) - Cvco2*(Pho[i] - Mu[i]*(1-Cstore[i])*E - pCstore[i]*E) #*X[i] not including X[i] 
        #since it is turbidostat and X is manitained constant (almost every hour dilution; )
        #see 07_OD_revised.XLSM
    else:
        dCstore = 0
        dCstore2 = - Cstoreloss 
        dX = 0
        dCO2 = Hco2*(CO2air - CO2[i])
    
    def upd(f,df):  #upd = update
        f[i+1] = f[i] + df*dt 
        return f
    
    Cstore = upd(Cstore,dCstore)
    Cstore2 = upd(Cstore2,dCstore2)
    X = upd(X,dX)
    CO2 = upd(CO2,dCO2)
    CO2 = CO2air*ones(size(o()))
    
    Qc[i+1] = QcBio + Cstore[i+1] 
    MuC[i] = (Qc[i+1]*X[i+1] - Qc[i]*X[i])/(Qc[i]*X[i]*dt)

#For the last step================
MuC[-1] = nan

MuC[t>=NightTime] = nan
CO2[t>=NightTime] = nan

#C:N ======================
#Average = 5.137
CN = 5.137*ones(size(Qc)) #C:N based on the data average
Ycn = CN/Qc #(mol mol-1) C:N of biomass
MuC0 = copy(MuC); MuC0[t>=NightTime] = 0
Vn = MuC0/CN #(molN molC-1 d-1) N uptake rate
CNdata = genfromtxt('..\\Data\\CNdataNO3.csv',delimiter=',').T

Vn[t>=NightTime] = nan
Ycn[t>=NightTime] = nan
CN[t>=NightTime ] = nan
Pho[t>=NightTime] = nan
#OOOOOOOOOOOOOOOOOOOOOOOOOOOO
# Data mpreparation
#OOOOOOOOOOOOOOOOOOOOOOOOOOOO

a = genfromtxt('C:\\Users\\keiin\\Google Drive\\Croco collaboration\\\
04 GAP\\02 Excel files\\Data\\DataNO3.csv',delimiter=',').T

OD = genfromtxt('C:\\Users\\keiin\\Google Drive\\Croco collaboration\\\
04 GAP\\02 Excel files\\Data\\ODdataNO3.csv',delimiter=',').T

QcData = a[24,2]    #(fg cell-1)

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# Plot
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

Xlabel = '$\mathit{t}$ (h)'

def sup(loc,Title,Ylabel):
    subplot(2,2,loc)
    title(Title,y=1.02)
    xlabel(Xlabel)
    ylabel(Ylabel)

def sup1(loc,Title,Ylabel):
    figure(loc)
   # title(Title,y=1.02)
    xlabel(Xlabel)
    ylabel(Ylabel)

def ebar(x,y,e):
    errorbar(x,y,(e,e),fmt='o',color='red',\
         elinewidth=1,markeredgecolor='k',ecolor='k',capthick=1,capsize=5)

def ebar2(x,y,e,l):
    errorbar(x,y,(e,e),fmt='o',color='red',\
         elinewidth=1,markeredgecolor='k',ecolor='k',capthick=1,capsize=5,label=str(l))

def Xlim():
    xlim(left=-0.3,right=14.3)
    xticks(arange(0,14.01,2))

def Xlim2():
    xlim(left=-0.3,right=24.3)
    xticks(arange(0,24.1,4))

def Xlim3():
    xlim(left=0,right=14)
    xticks(arange(0,14.01,2))

def sf(name):
    Savefig3('02\\02 GAP1\\04 paper\\NO3',name,300)


# #figure(1,figsize=(12,9))
# 
sup1(1,'$\mathit{C_{Sto}}$','$\mathit{C_{Sto}}$ (mol C mol C$^{-1}$)')
plot(th,Cstore)
plot(th,Cstore2,'b--')
ebar(a[12],a[13]/QcData,a[14]/QcData)
Xlim2()
ylim(0.,1.1)
sf('Cstore')
 
sup1(2,'$\mathit{DIC}$','[$\mathit{DIC}]$ ($\mu$M C)')
plot(th, CO2)
ebar(a[4],a[5],a[6])
ylim()
Xlim()
ylim(-100,1000)
sf('CO2')
 
sup1(3,'$\mathit{XQ_C}$','[$\mathit{C_{Cell}}$]')
plot(th,X*Qc/(X[0]*(Cstore[0]+1)))
plot(OD[0],OD[1])
Xlim2()
ylim(0.97,1.5)
sf('OD')
 
sup1(4,'$\mathit{\mu_{C}}$','$\mathit{\mu_{C}}$ (d$^{-1}$)')
plot(th,MuC)
ebar(a[0],a[1]*24,a[2]*24)
ylim(-0.2,3.0)
Xlim()
sf('MuC')
 
rcParams.update({'xtick.top': 'False'})
rcParams.update({'ytick.direction': 'in'})
rcParams.update({'xtick.direction': 'out'})
sup1(5,'','Fate of C (d$^{-1}$)')
stackplot(th,Mu*1,pCstore,(pCstore-Mu*Cstore)*E+Mu*E)
Xlim3()
ylim(0,5)

Names = ["Biomass","C store","Biosyn. cost"]
Colors = ['blue','red','green']
stackplot([],[],[],[],colors=Colors[::-1],labels=Names[::-1])
legend(loc=4,edgecolor = 'k')
sf('C_fate')
# 
# figure(6)
# plot(th,pCstore)
# 
# figure(7)
# plot(th,Cstore)

sup1(8,'','C:N (mol mol$^{-1}$)')
plot(th,CN,label='Model')
ebar2(CNdata[0],CNdata[1],CNdata[2],"Data")
#plot(CNdata[0],CNdata[1],'o',label='Data')
#plot(th,Ycn,label='Bio')
legend(loc=4,edgecolor='k')
ylim(0,8)
sf('CN')

sup1(9,'','V$_{N}$ (molN molC$^{-1}$ d$^{-1}$)')
plot(th,Vn*(QcBio + Cstore))
ylim(bottom=0)
print(Vn)
sf('Vn')

sup1(10,'','e$_{Cfix}$ : e$_{NO3}$ electron ratio')
#plot(th,CN*4/8)
plot(th,MuC0*4/(Vn*8))  
print(MuC0*4/(Vn*8)) 
ylim(0,5)
sf('CNelectron')

sup1(11,'','Cell N (molN molC$^{-1}$)')
#plot(th,CN,label='Model')
plot(th,1/Ycn)#,label='Bio')
#legend(loc=4,edgecolor='k')
ylim(0,0.4)
sf('N_Cell')

#ylim(0,8)

#Savefig3('02\\02 GAP1\\02 testing','All',300)
#Savetxt(Mu,'02\\02 GAP1\\02 testing','DIC')

show()





