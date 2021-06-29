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
E = Enstore
#E = 0

#----------------------
# Storage
#----------------------
pCstoreMax = 8 #(d-1)
CstoreMax = 0.92  #(molC molC-1)

#:::::::::::::::::::::::::::::::::::::::::::::::
# Key parameters; These are really keys
#:::::::::::::::::::::::::::::::::::::::::::::::
PhoMax = 4.5    #(d-1) Maximum photosynthetsis rate
Kco2 = 100     #(uM)

Hco2 = 2 #(d-1) gas exchange constant 
Cvco2 = 2.9e3 #(dimensionless) factor for adjusting CO2 consumption

Cn2Fix = 1.3 #(d-1)

Cstoreloss = 1.7 #(d-1) Cstoreloss due to its conversion to N containing molecules (biomass)

#============================
# 
#============================
dts = 10  #(s) dt in second
dt = dts/86400 #(d) dt in day
t = arange(0,86400/24*24+dts,dts)/86400 #(d) time
th = t*24 #(h) time
U = arange(size(t))
NightTime = 3600*14/86400
InitN2fix = 3600*15.5/86400
EndN2fix = 3600*19.5/86400
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
CO2[0] = 672 #(uM)
CO2 = CO2air*ones(size(o()))
Cstore[0] = 0.38 #(molC molC-1)
Cstore2[0] = Cstore[0]
X[0] = 1

#=====================

for i in U[:-1]:
    if t[i] < NightTime:
        Pho[i] = PhoMax*CO2[i]/(CO2[i]+Kco2)
        pCstore[i] = pCstoreMax*(CstoreMax - Cstore[i])/CstoreMax
        Mu[i] = (Pho[i] - pCstore[i])/(1+E)
    
    if Pho[i] <= 0:
        Pho[i] = 0
        pCstore[i] = 0
        Mu[i] = 0
    
    if Pho[i] < pCstore[i]:
        pCstore[i] = Pho[i]
        Mu[i] = 0
    
    Qc[i] = QcBio + Cstore[i]
    
    if t[i]<NightTime:
    
        dCstore = Pho[i] - Mu[i]*Qc[i]*(1+E)
        dCstore = pCstore[i] - Mu[i]*Cstore[i]
        dCstore2 = dCstore
        
        dX = Mu[i]*X[i]#/Qc[i]
        dCO2 = Hco2*(CO2air - CO2[i]) - Cvco2*(Pho[i] - Mu[i]*E) #*X[i] not including X[i] 
        #since it is turbidostat and X is manitained constant (almost every hour dilution; )
        #see 07_OD_revised.XLSM
    
    elif t[i]>=InitN2fix and t[i]<EndN2fix:
        dCstore = - Cn2Fix
        dCstore2 = - Cn2Fix - Cstoreloss
        Mu[i] = 0
    
    else:
        dCstore = 0
        dCstore2 = 0
        dX = 0
        dCO2 = 0
        Mu[i] = 0
    
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
    MuC[i] = (pCstore[i]+Mu[i])/Qc[i]

   # MuC[i] = (Pho[i] - Mu[i]*E)/Qc[i]

#For the last step================
MuC[-1] = nan

MuC[t>=NightTime] = nan
CO2[t>=NightTime] = nan
pCstore[t>=NightTime] = nan

#C:N ======================
Ycn = 3.86 #(mol mol-1) C:N of biomass
CN = Qc*Ycn
Vn = Mu/Ycn #(molN molC d-1)
iNight = int(NightTime/dt)
SigmaVn = (X[iNight] - X[0])/Ycn #(molN molC d-1)
N2fix = SigmaVn/(EndN2fix-InitN2fix)*ones(size(Vn)) #(molN molC d-1)
N2fix[t<InitN2fix] = 0
N2fix[t>EndN2fix] = 0
CNdata = genfromtxt('..\\Data\\CNdataN2.csv',delimiter=',').T
CN[t>=NightTime] = nan
#OOOOOOOOOOOOOOOOOOOOOOOOOOOO
# Data preparation
#OOOOOOOOOOOOOOOOOOOOOOOOOOOO

a = genfromtxt('C:\\Users\\keiin\\Google Drive\\Croco collaboration\\\
04 GAP\\02 Excel files\\Data\\DataN2.csv',delimiter=',').T

OD = genfromtxt('C:\\Users\\keiin\\Google Drive\\Croco collaboration\\\
04 GAP\\02 Excel files\\Data\\ODdataN2.csv',delimiter=',').T

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
    Savefig3('02\\02 GAP1\\04 paper\\N2fix',name,300)

#figure(1,figsize=(12,9))

sup1(1,'$\mathit{C_{Store}}$','$\mathit{C_{Sto}}$ (mol C mol C$^{-1}$)')
plot(th,Cstore,label='Model')
plot(th,Cstore2,'b--')
ebar2(a[12],a[13]/QcData,a[14]/QcData,'Data')
ylim(bottom=0)
Xlim2()
ylim(0.,1.1)
legend(loc=4,edgecolor='k')
sf('Cstore')
#Xlim()

sup1(2,'$\mathit{DIC}$','[$\mathit{DIC}$] ($\mu$M C)')
plot(th, CO2)
ebar(a[4],a[5],a[6])
Xlim()
ylim(-100,1000)
sf('CO2')

sup1(3,'$\mathit{XQ_C}$','[$\mathit{C_{Cell}}$]')
plot(th,X*Qc/(X[0]*(Cstore[0]+1)),label='Model')
plot(OD[0],OD[1],label='Data')
Xlim2()
ylim(0.97,1.5)
legend(loc=4,edgecolor='k')
sf('OD')

sup1(4,'$\mathit{\mu Q_C}$','$\mathit{\mu_{C}}$ (d$^{-1}$)')
plot(th,MuC,label='Model')
ebar2(a[0],a[1]*24,a[2]*24,'Data')
ylim(bottom = -0.2)
Xlim()
ylim(-0.2,3.0)
legend(loc=1,edgecolor='k')
sf('MuC')

rcParams.update({'xtick.top': 'False'})
rcParams.update({'ytick.direction': 'in'})
rcParams.update({'xtick.direction': 'out'})
sup1(5,'','Fate of C (d$^{-1}$)')
Names = ["Biomass","C store","Biosyn. cost"]
Colors = ['blue','red','green']
stackplot([],[],[],[],colors=Colors[::-1],labels=Names[::-1])
stackplot(th,Mu,pCstore,Mu*E,colors=Colors)
Xlim3()
ylim(0,5)
legend(loc=4,edgecolor = 'k')
sf('C_fate')

sup1(6,'','C:N (mol mol$^{-1}$)')
plot(th,CN,label='Model')
#plot(CNdata[0],CNdata[1],'o',label='Data')
ebar2(CNdata[0],CNdata[1],CNdata[2],'Data')
#plot(th,Ycn*CN/CN,label='Bio',zorder=0)
legend(loc=4,edgecolor='k')
ylim(0,8)
sf('CN')


sup1(7,'','N$_{2}$ fix (molN molC$^{-1}$)')
plot(th,N2fix)
#Savefig3('02\\02 GAP1\\02 testing','All',300)
#Savetxt(Mu,'02\\02 GAP1\\02 testing','DIC')

sup1(11,'','Cell N (molN molC$^{-1}$)')
#plot(th,CN,label='Model')
plot(th,1/Ycn*CN/CN)#,label='Bio')
#legend(loc=4,edgecolor='k')
ylim(0,0.4)
sf('CN_Bio')

show()





