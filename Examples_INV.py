# -*- coding: utf-8 -*-
"""
Created on Sun May 12 19:59:48 2019

@author: ollenil
"""

# Commentare il codice in italiano, così da distinguere i nostri commenti dai (pochi) commenti del creatore :D

import numpy as np
import math
import csv
from future.backports.http.cookiejar import iso2time
from Analysis_Method import LOAD, PANEL, ASOLVE
from Inverse_Method import ISOLVE
from Post_noLatex_INV import CP, PLOT_CP
from Post_noLatex_INV import FIELD_U_V_CP, STREAMLINES
from Post_noLatex_INV import FIELD_PLOT_CP_U_V,FIELD_PLOT_STREAMLINES

# Dati test Designers
"""
from Aerofoil_Generation import NACA4series
thickness=0.12
M=120
c=1
alpha=15
GE=False
h=0.0
save=True
savename = 'NACA0012Example'
"""

# Test Designers
"""
def Analyse_NACA(thickness,M,c,alpha,GE,h):
    XP,ZP=NACA4series(thickness,M=120,c=1,plot=True)
    GAMMA,PSI,CL=ASOLVE(XP,ZP,alpha,UINF=1,GE=GE,h=0.2)
    CPX=CP(GAMMA,XP)
    PLOT_CP(CPX,XP,ZP,GE=GE,save=save,savename=savename)
"""    
    
def Analyse_txtfile(filepath,filename,M,alpha,GE,h):
    XC,ZC=LOAD(filepath,filename,plot=False)
    XP,ZP=PANEL(XC,ZC,M=M,plot=False)
    GAMMA,PSI,CL=ASOLVE(XP,ZP,alpha,UINF=1,GE=GE,h=0.2)
    CPX=CP(GAMMA,XP)
    return XP,ZP,GAMMA,PSI,CPX

def Plot_Field(xlimit,zlimit,NX,NZ,filepath,filename,M,alpha,GE,h):
    XC,ZC=LOAD(filepath,filename,plot=False)
    XP,ZP=PANEL(XC,ZC,M=M,plot=False)
    GAMMA,PSI,CL=ASOLVE(XP,ZP,alpha,UINF=1,GE=GE,h=0.2)
    FIELD_PSI=np.array((0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9,1.0))
    PSIFIELD=np.concatenate((np.ravel(PSI),FIELD_PSI))
    
    #Calculate velocity and Cp field
    (U, W, VM, CPF, XG, ZG)=FIELD_U_V_CP(xlimit,zlimit,NX,NZ,GAMMA,XP,ZP,alpha,h)
    #plot velocity and Cp field
    FIELD_PLOT_CP_U_V(U, W, VM, CPF, XG, ZG, XP, ZP,zlimit, h=h,mirror=False)
    #Calculate Stramlines,    
    XZPSI=STREAMLINES(xlimit,zlimit,NX,NZ,GAMMA,XP,ZP,alpha,PSIFIELD,h)
    #plot 
    FIELD_PLOT_STREAMLINES(XZPSI, XP, ZP,zlimit, h=h,mirror=False,savename='test',save=True)
    #plot(CL)

filepath='D:\\SC_Polito\\SC_23\\Progetti\\External_CFD\\Panels_Method\\' # Inserisci la tua path in questo formato
filename='Inverted_NRL_071_coordinates.txt'   
M=120
alpha=0
UINF = 1
GE=True
h=0.2
save=True
savename = 'inverse_test'

# Griglia di calcolo attorno ai profili per plottare i campi di moto            
xlimit=[-0.5,1.8]
zlimit=[-0.1,1.0]
NX=40
NZ=40   

# Calcolo di XP, CP, GAMMA, ..., da utilizzare poi come input per il metodo inverso
XP, ZP, GAMMAI, PSII, CPX = Analyse_txtfile(filepath,filename,M,alpha,GE,h)

# Questa parte di codice serve a generare un file CSV con colonne XP, CP del profilo originale dato in input
# Può essere utile se si vuole modificare manualmente la distribuzione di CP dal file CSV
"""
with open('CP.csv','w',newline='') as file:
    writer = csv.writer(file)
    for row in range(len(XP[0])):
        writer.writerow([f"{float(XP[0][row])}", f"{float(CPX[0][row])}"])
"""

# Questa parte di codice legge in input un file CSV con colonne XP, CP obiettivo, cioè quello che si vuole in output
# In questo primo test è stata modificata solo la distribuzione di CP del main. Non sarà comunque complicato
# generalizzare il codice a tutti i casi.
XSP = XP
CP_SP = [[],[]]
with open('CP_NEW.csv','r') as file:
    reader = csv.reader(file)
    for row in reader:
        CP_SP[0] = np.append(CP_SP[0], float(row[1])) # Modifica del CP sul main
CP_SP[1] = CPX[1] # Distribuzione di CP sul flap uguale a quella originale

# GAMMA obiettivo ricavato dal CP obiettivo tramite l'equazione 22 del paper
GAMMASP = [[],[],[],[]]
for i in CP_SP[0]:
    GAMMASP[0] = np.append(GAMMASP[0], float(UINF * math.sqrt(1-i)))
GAMMASP[1] = GAMMAI[1]
K = 2 # Generalizzare per K (numero di profili) qualsiasi
for k in range(K):
    for i in GAMMASP[k]:
          GAMMASP[k+K] = np.append(GAMMASP[k+K], i)

# Limiti dei segmenti di design
iS = [1, 11] 
iE = [10, 15]

# Calcolo inverso
isolve_sol = []
isolve_sol = ISOLVE(XP,ZP,GAMMAI,PSII,GAMMASP,XSP,iS,iE,ALPHA=0,UINF=1,GE=True,h=0.2)
XNEW = isolve_sol[0]
ZNEW = isolve_sol[1]
GAMMANEW = isolve_sol[2]
PSISOL = isolve_sol[3]
CPXNEW = CP(GAMMANEW,XNEW)
PLOT_CP(CPX,XP,ZP,CPXNEW,XNEW,ZNEW,GE=GE,save=save,savename=savename)