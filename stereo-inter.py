#!/usr/bin/python
from __future__ import division
import numpy as np
from Tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib import pyplot as plt
from PIL import Image
import PngImagePlugin
import ttk
import sys
import os
from fractions import Fraction
from tkFileDialog import *
#import pdb

pi=3.14159265359


#Image._initialized=2
###################################################################"
##### Fonction projection  sur l'abaque
####################################################################

def proj(x,y,z): 
  
    if z==1: 
        X=0
        Y=0
    elif z<-0.000001:
        X=250
        Y=250
    else: 
            
        X=x/(1+z)
        Y=y/(1+z)
    
    return np.array([X,Y],float) 
    
###################################################################"
##### Fonction rotation 
####################################################################

def rotation(phi1,phi,phi2):
   phi1=phi1*pi/180;
   phi=phi*pi/180;
   phi2=phi2*pi/180;
   R=np.array([[np.cos(phi1)*np.cos(phi2)-np.cos(phi)*np.sin(phi1)*np.sin(phi2),
            -np.cos(phi)*np.cos(phi2)*np.sin(phi1)-np.cos(phi1)*
            np.sin(phi2),np.sin(phi)*np.sin(phi1)],[np.cos(phi2)*np.sin(phi1)
            +np.cos(phi)*np.cos(phi1)*np.sin(phi2),np.cos(phi)*np.cos(phi1)
            *np.cos(phi2)-np.sin(phi1)*np.sin(phi2), -np.cos(phi1)*np.sin(phi)],
            [np.sin(phi)*np.sin(phi2), np.cos(phi2)*np.sin(phi), np.cos(phi)]],float)
   return R

####################################################################
##### Fonction rotation autour d'un axe 
####################################################################

def Rot(th,a,b,c):
   th=th*pi/180;
   aa=a/np.linalg.norm([a,b,c]);
   bb=b/np.linalg.norm([a,b,c]);
   cc=c/np.linalg.norm([a,b,c]);
   c1=np.array([[1,0,0],[0,1,0],[0,0,1]],float)
   c2=np.array([[aa**2,aa*bb,aa*cc],[bb*aa,bb**2,bb*cc],[cc*aa,
                cc*bb,cc**2]],float)
   c3=np.array([[0,-cc,bb],[cc,0,-aa],[-bb,aa,0]],float)
   R=np.cos(th)*c1+(1-np.cos(th))*c2+np.sin(th)*c3

   return R    


####################################################################
##### Fonction cristal
####################################################################
def cristA():
    global axesA,axeshA,DA,DstarA,VA
    a=eval(aA_entry.get())
    b=eval(bA_entry.get())
    c=eval(cA_entry.get())
    alp=eval(alpA_entry.get())
    bet=eval(betA_entry.get())
    gam=eval(gamA_entry.get())
    e=eval(eA_entry.get())
    d2=eval(dA_label_var.get())
    alp=alp*pi/180;
    bet=bet*pi/180;
    gam=gam*pi/180;
    VA=a*b*c*np.sqrt(1-(np.cos(alp)**2)-(np.cos(bet))**2-(np.cos(gam))**2+2*b*c*np.cos(alp)*np.cos(bet)*np.cos(gam))
    DA=np.array([[a,b*np.cos(gam),c*np.cos(bet)],[0,b*np.sin(gam),  c*(np.cos(alp)-np.cos(bet)*np.cos(gam))/np.sin(gam)],[0,0,VA/(a*b*np.sin(gam))]])
    DstarA=np.transpose(np.linalg.inv(DA))
    G=np.array([[a**2,a*b*np.cos(gam),a*c*np.cos(bet)],[a*b*np.cos(gam),b**2,b*c*np.cos(alp)],[a*c*np.cos(bet),b*c*np.cos(alp),c**2]])    
    axes=np.zeros(((2*e+1)**3-1,3))
    axesh=np.zeros(((2*e+1)**3-1,4))
    id=0
    for i in range(-e,e+1):
        for j in range(-e,e+1):
            for k in range(-e,e+1):
                if (i,j,k)!=(0,0,0):
                    d=1/(np.sqrt(np.dot(np.array([i,j,k]),np.dot(np.linalg.inv(G),np.array([i,j,k])))))
                    if d>d2*0.1*np.amax([a,b,c]):
                        if var_uvw.get()==0:                    
                            Ma=np.dot(DstarA,np.array([i,j,k],float))
                            axesh[id,0]=Ma[0]
                            axesh[id,1]=Ma[1]
                            axesh[id,2]=Ma[2]
                            axesh[id,3]=0
                            axes[id,:]=np.array([i,j,k],float)
                        else:
                            Ma=np.dot(DA,np.array([i,j,k],float))
                            axesh[id,0]=Ma[0]
                            axesh[id,1]=Ma[1]
                            axesh[id,2]=Ma[2]
                            axesh[id,3]=1
                            axes[id,:]=np.array([i,j,k],float)
                        id=id+1
    axesA=axes
    
    axeshA=axesh
                 
    return axesA,axeshA,DA,DstarA,VA

def cristB():
    global axesB,axeshB,DB,DstarB,VB
    a=eval(aB_entry.get())
    b=eval(bB_entry.get())
    c=eval(cB_entry.get())
    alp=eval(alpB_entry.get())
    bet=eval(betB_entry.get())
    gam=eval(gamB_entry.get())
    e=eval(eB_entry.get())
    d2=eval(dB_label_var.get())
    alp=alp*pi/180;
    bet=bet*pi/180;
    gam=gam*pi/180;
    VB=a*b*c*np.sqrt(1-(np.cos(alp)**2)-(np.cos(bet))**2-(np.cos(gam))**2+2*b*c*np.cos(alp)*np.cos(bet)*np.cos(gam))
    DB=np.array([[a,b*np.cos(gam),c*np.cos(bet)],[0,b*np.sin(gam),  c*(np.cos(alp)-np.cos(bet)*np.cos(gam))/np.sin(gam)],[0,0,VB/(a*b*np.sin(gam))]])
    DstarB=np.transpose(np.linalg.inv(DB))
    G=np.array([[a**2,a*b*np.cos(gam),a*c*np.cos(bet)],[a*b*np.cos(gam),b**2,b*c*np.cos(alp)],[a*c*np.cos(bet),b*c*np.cos(alp),c**2]])    
    axes=np.zeros(((2*e+1)**3-1,3))
    axesh=np.zeros(((2*e+1)**3-1,4))
    id=0
    for i in range(-e,e+1):
        for j in range(-e,e+1):
            for k in range(-e,e+1):
                if (i,j,k)!=(0,0,0):
                    d=1/(np.sqrt(np.dot(np.array([i,j,k]),np.dot(np.linalg.inv(G),np.array([i,j,k])))))
                    if d>d2*0.1*np.amax([a,b,c]):
                             if var_uvw.get()==0:                    
                                Ma=np.dot(DstarB,np.array([i,j,k],float))
                                axesh[id,0]=Ma[0]
                                axesh[id,1]=Ma[1]
                                axesh[id,2]=Ma[2]
                                axesh[id,3]=0
                                axes[id,:]=np.array([i,j,k],float)
                             else:
                                Ma=np.dot(DB,np.array([i,j,k],float))
                                axesh[id,0]=Ma[0]
                                axesh[id,1]=Ma[1]
                                axesh[id,2]=Ma[2]
                                axesh[id,3]=1
                                axes[id,:]=np.array([i,j,k],float)
                             id=id+1
    axesB=axes
    
    axeshB=axesh
                 
    return axesB,axeshB,DB,DstarB,VB
    
    
def dmA():
    global dmipA
    a=f.add_subplot(111)    
    a.figure.clear()
    a=f.add_subplot(111)      
    dmipA=dmipA-eval(dA_entry.get())
    dA_label_var.set(dmipA)
    cristA()
    trace()
        
    return dmipA
    
def dmB():
    global dmipB
    a=f.add_subplot(111)    
    a.figure.clear()
    a=f.add_subplot(111)      
    dmipB=dmipB-eval(dB_entry.get())
    dB_label_var.set(dmipB)
    cristB()
    trace()
    
    return dmipB

def dpA():
    global dmipA
    a=f.add_subplot(111)    
    a.figure.clear()
    a=f.add_subplot(111)      
    dmipA=dmipA+eval(dA_entry.get())
    dA_label_var.set(dmipA)    
    cristA()    
    trace()
       
    return dmipA 

def dpB():
    global dmipB
    a=f.add_subplot(111)    
    a.figure.clear()
    a=f.add_subplot(111)      
    dmipB=dmipB+eval(dB_entry.get())
    dB_label_var.set(dmipB)    
    cristB()    
    trace()
       
    return dmipB 
    

####################################################################
##### Fonction ajouter un pole
####################################################################
def poleA(pole1,pole2,pole3):
    global MA,axesA,axeshA,Ta,VA,DA,DstarA
    
        
    if var_hexaA.get()==1:
        if var_uvw.get()==1:
            pole1a=2*pole1+pole2
            pole2a=2*pole2+pole1
            pole1=pole1a
            pole2=pole2a
    
    Gs=np.array([pole1,pole2,pole3],float)
    #print(Gs)
    if var_uvw.get()==0:                    
            Gsh=np.dot(DstarA,Gs)/np.linalg.norm(np.dot(DstarA,Gs))
    else:
        Gsh=np.dot(DA,Gs)/np.linalg.norm(np.dot(DA,Gs))
     
    S=np.dot(MA,Gsh)
    if S[2]<0:
        S=-S
        Gsh=-Gsh
        pole1=-pole1
        pole2=-pole2
        pole3=-pole3
   
    
    axesA=np.vstack((axesA,np.array([pole1,pole2,pole3])))
    axesA=np.vstack((axesA,np.array([-pole1,-pole2,-pole3])))
    Ta=np.vstack((Ta,np.array([S[0],S[1],S[2]])))
    Ta=np.vstack((Ta,np.array([-S[0],-S[1],-S[2]])))
    axeshA=np.vstack((axeshA,np.array([Gsh[0],Gsh[1],Gsh[2],0])))
    axeshA=np.vstack((axeshA,np.array([-Gsh[0],-Gsh[1],-Gsh[2],1])))
    return axesA,axeshA,Ta
    
def poleB(pole1,pole2,pole3):
    global MB,axesB,axeshB,Tb,VB,DB,DstarB
    
    
    
    if var_hexaB.get()==1:
        if var_uvw.get()==1:
            pole1a=2*pole1+pole2
            pole2a=2*pole2+pole1
            pole1=pole1a
            pole2=pole2a
    
    Gs=np.array([pole1,pole2,pole3],float)
    #print(Gs)
    if var_uvw.get()==0:                    
            Gsh=np.dot(DstarB,Gs)/np.linalg.norm(np.dot(DstarB,Gs))
    else:
        Gsh=np.dot(DB,Gs)/np.linalg.norm(np.dot(DB,Gs))
     
    S=np.dot(MB,Gsh)
    if S[2]<0:
        S=-S
        Gsh=-Gsh
        pole1=-pole1
        pole2=-pole2
        pole3=-pole3
     
    
    axesB=np.vstack((axesB,np.array([pole1,pole2,pole3])))
    axesB=np.vstack((axesB,np.array([-pole1,-pole2,-pole3])))
    Tb=np.vstack((Tb,np.array([S[0],S[1],S[2]])))
    Tb=np.vstack((Tb,np.array([-S[0],-S[1],-S[2]])))
    axeshB=np.vstack((axeshB,np.array([Gsh[0],Gsh[1],Gsh[2],0])))
    axeshB=np.vstack((axeshB,np.array([-Gsh[0],-Gsh[1],-Gsh[2],1])))
    return axesB,axeshB,Tb
    
    
def addpoleA_sym():
        
    pole1A=eval(pole1A_entry.get())
    pole2A=eval(pole2A_entry.get())
    pole3A=eval(pole3A_entry.get())
    poleA(pole1A,pole2A,pole3A)
    poleA(pole1A,pole2A,-pole3A)
    poleA(pole1A,-pole2A,pole3A)
    poleA(-pole1A,pole2A,pole3A)
    poleA(pole2A,pole1A,pole3A)
    poleA(pole2A,pole1A,-pole3A)
    poleA(pole2A,-pole1A,pole3A)
    poleA(-pole2A,pole1A,pole3A)
    poleA(pole2A,pole3A,pole1A)
    poleA(pole2A,pole3A,-pole1A)
    poleA(pole2A,-pole3A,pole1A)
    poleA(-pole2A,pole3A,pole1A)
    poleA(pole1A,pole3A,pole2A)
    poleA(pole1A,pole3A,-pole2A)
    poleA(pole1A,-pole3A,pole2A)
    poleA(-pole1A,pole3A,pole2A)
    poleA(pole3A,pole1A,pole2A)
    poleA(pole3A,pole1A,-pole2A)
    poleA(pole3A,-pole1A,pole2A)
    poleA(-pole3A,pole1A,pole2A)
    poleA(pole3A,pole2A,pole1A)
    poleA(pole3A,pole2A,-pole1A)
    poleA(pole3A,-pole2A,pole1A)
    poleA(-pole3A,pole2A,pole1A)
    trace()

def addpoleB_sym():
        
    pole1B=eval(pole1B_entry.get())
    pole2B=eval(pole2B_entry.get())
    pole3B=eval(pole3B_entry.get())
    poleB(pole1B,pole2B,pole3B)
    poleB(pole1B,pole2B,-pole3B)
    poleB(pole1B,-pole2B,pole3B)
    poleB(-pole1B,pole2B,pole3B)
    poleB(pole2B,pole1B,pole3B)
    poleB(pole2B,pole1B,-pole3B)
    poleB(pole2B,-pole1B,pole3B)
    poleB(-pole2B,pole1B,pole3B)
    poleB(pole2B,pole3B,pole1B)
    poleB(pole2B,pole3B,-pole1B)
    poleB(pole2B,-pole3B,pole1B)
    poleB(-pole2B,pole3B,pole1B)
    poleB(pole1B,pole3B,pole2B)
    poleB(pole1B,pole3B,-pole2B)
    poleB(pole1B,-pole3B,pole2B)
    poleB(-pole1B,pole3B,pole2B)
    poleB(pole3B,pole1B,pole2B)
    poleB(pole3B,pole1B,-pole2B)
    poleB(pole3B,-pole1B,pole2B)
    poleB(-pole3B,pole1B,pole2B)
    poleB(pole3B,pole2B,pole1B)
    poleB(pole3B,pole2B,-pole1B)
    poleB(pole3B,-pole2B,pole1B)
    poleB(-pole3B,pole2B,pole1B)
    trace()
        
def addpoleA():
    
    pole1A=eval(pole1A_entry.get())
    pole2A=eval(pole2A_entry.get())
    pole3A=eval(pole3A_entry.get())
    poleA(pole1A,pole2A,pole3A)
    trace()
def addpoleB():
    
    pole1B=eval(pole1B_entry.get())
    pole2B=eval(pole2B_entry.get())
    pole3B=eval(pole3B_entry.get())
    poleB(pole1B,pole2B,pole3B)
    trace()
####################################################################
##### Fonction tracer plan
####################################################################
def trace_planA(pole1,pole2,pole3):
    global MA,axesA,axeshA,T,VA,DA,DstarA,trPA
    
    
    pole_i=0
    
    
    if var_hexaA.get()==1:
        if var_uvw.get()==1:
            pole1=2*eval(pole1A_entry.get())+eval(pole2A_entry.get())
            pole2=2*eval(pole2A_entry.get())+eval(pole1A_entry.get())
            pole3=eval(pole3A_entry.get())
            pole_i=1

        
    trPA=np.vstack((trPA,np.array([pole1,pole2,pole3,pole_i])))
    b=np.ascontiguousarray(trPA).view(np.dtype((np.void, trPA.dtype.itemsize * trPA.shape[1])))
    
    trPA=np.unique(b).view(trPA.dtype).reshape(-1, trPA.shape[1])
    #print(trP)
    
def trace_planB(pole1,pole2,pole3):
    global MB,axesB,axeshB,T,VB,DB,DstarB,trPB
    
    
    pole_i=0
    
    
    if var_hexaB.get()==1:
        if var_uvw.get()==1:
            pole1=2*eval(pole1B_entry.get())+eval(pole2B_entry.get())
            pole2=2*eval(pole2B_entry.get())+eval(pole1B_entry.get())
            pole3=eval(pole3B_entry.get())
            pole_i=1

        
    trPB=np.vstack((trPB,np.array([pole1,pole2,pole3,pole_i])))
    b=np.ascontiguousarray(trPB).view(np.dtype((np.void, trPB.dtype.itemsize * trPB.shape[1])))
    
    trPB=np.unique(b).view(trPB.dtype).reshape(-1, trPB.shape[1])
    #print(trP)    
    
def trace_addplanA():
    global MA,axesA,axeshA,T,VA,DA,DstarA,trPA
    
    pole1=eval(pole1A_entry.get())
    pole2=eval(pole2A_entry.get())
    pole3=eval(pole3A_entry.get())
    
    trace_planA(pole1,pole2,pole3)
    trace()

def trace_addplanB():
    global MB,axesB,axeshB,T,VB,DB,DstarB,trPB
    
    pole1=eval(pole1B_entry.get())
    pole2=eval(pole2B_entry.get())
    pole3=eval(pole3B_entry.get())
    
    trace_planB(pole1,pole2,pole3)
    trace()    

def trace_plan_symA():
    global MA,axesA,axeshA,T,VA,DA,DstarA,GA    
    pole1=eval(pole1_entry.get())
    pole2=eval(pole2_entry.get())
    pole3=eval(pole3_entry.get())
    a=eval(aA_entry.get())
    b=eval(bA_entry.get())
    c=eval(cA_entry.get())
    alp=eval(alpA_entry.get())
    bet=eval(betA_entry.get())
    gam=eval(gamA_entry.get())
    alp=alp*pi/180;
    bet=bet*pi/180;
    gam=gam*pi/180;
    GA=np.array([[a**2,a*b*np.cos(gam),a*c*np.cos(bet)],[a*b*np.cos(gam),b**2,b*c*np.cos(alp)],[a*c*np.cos(bet),b*c*np.cos(alp),c**2]])     
    v=d(pole1,pole2,pole3)
    
    trace_planA(pole1,pole2,pole3)
    if np.abs(alp-pi/2)<0.001 and np.abs(bet-pi/2)<0.001 and np.abs(gam-2*pi/3)<0.001:
        trace_plan(pole1,pole2,pole3)
        trace_plan(pole1,pole2,-pole3)
        trace_plan(pole2,pole1,pole3)
        trace_plan(pole2,pole1,-pole3)
        trace_plan(-pole1-pole2,pole2,pole3)
        trace_plan(-pole1-pole2,pole2,-pole3)
        trace_plan(pole1,-pole1-pole2,pole3)
        trace_plan(pole1,-pole1-pole2,-pole3)
        trace_plan(pole2,-pole1-pole2,pole3)
        trace_plan(pole2,-pole1-pole2,-pole3)
        trace_plan(-pole1-pole2,pole1,pole3)
        trace_plan(-pole1-pole2,pole1,-pole3)

    else:
        if np.abs(d(pole1,pole2,-pole3)-v)<0.001:
                trace_plan(pole1,pole2,-pole3)
        if np.abs(d(pole1,-pole2,pole3)-v)<0.001:
                trace_plan(pole1,-pole2,pole3)
        if np.abs(d(-pole1,pole2,pole3)-v)<0.001:
            trace_plan(-pole1,pole2,pole3)
        if np.abs(d(pole2,pole1,pole3)-v)<0.001:
            trace_plan(pole2,pole1,pole3)
        if np.abs(d(pole2,pole1,-pole3)-v)<0.001:
            trace_plan(pole2,pole1,-pole3)
        if np.abs(d(pole2,-pole1,pole3)-v)<0.001:
            trace_plan(pole2,-pole1,pole3)
        if np.abs(d(-pole2,pole1,pole3)-v)<0.001:
            trace_plan(-pole2,pole1,pole3)
        if np.abs(d(pole2,pole3,pole1)-v)<0.001:
            trace_plan(pole2,pole3,pole1)
        if np.abs(d(pole2,pole3,pole1)-v)<0.001:
            trace_plan(pole2,pole3,-pole1)
        if np.abs(d(pole2,-pole3,pole1)-v)<0.001:
            trace_plan(pole2,-pole3,pole1)
        if np.abs(d(-pole2,pole3,pole1)-v)<0.001:
            trace_plan(-pole2,pole3,pole1)
        if np.abs(d(pole1,pole3,pole2)-v)<0.001:
            trace_plan(pole1,pole3,pole2)
        if np.abs(d(pole1,pole3,-pole2)-v)<0.001:
            trace_plan(pole1,pole3,-pole2)
        if np.abs(d(pole1,-pole3,pole2)-v)<0.001:
            trace_plan(pole1,-pole3,pole2)
        if np.abs(d(-pole1,pole3,pole2)-v)<0.001:
            trace_plan(-pole1,pole3,pole2)
        if np.abs(d(pole3,pole1,pole2)-v)<0.001:
            trace_plan(pole3,pole1,pole2)
        if np.abs(d(pole3,pole1,-pole2)-v)<0.001:
            trace_plan(pole3,pole1,-pole2)
        if np.abs(d(pole3,-pole1,pole2)-v)<0.001:
            trace_plan(pole3,-pole1,pole2)
        if np.abs(d(-pole3,pole1,pole2)-v)<0.001:
            trace_plan(-pole3,pole1,pole2)
        if np.abs(d(pole3,pole2,pole1)-v)<0.001:
            trace_plan(pole3,pole2,pole1)
        if np.abs(d(pole3,pole2,-pole1)-v)<0.001:
            trace_plan(pole3,pole2,-pole1)
        if np.abs(d(pole3,-pole2,pole1)-v)<0.001:
            trace_plan(pole3,-pole2,pole1)
        if np.abs(d(pole3,pole2,pole1)-v)<0.001:
            trace_plan(pole3,pole2,pole1)
    trace()
def trace_plan_symB():
    global MB,axesB,axeshB,T,VB,DB,DstarB,GB    
    pole1=eval(pole1_entry.get())
    pole2=eval(pole2_entry.get())
    pole3=eval(pole3_entry.get())
    a=eval(aB_entry.get())
    b=eval(bB_entry.get())
    c=eval(cB_entry.get())
    alp=eval(alpB_entry.get())
    bet=eval(betB_entry.get())
    gam=eval(gamB_entry.get())
    alp=alp*pi/180;
    bet=bet*pi/180;
    gam=gam*pi/180;
    GB=np.array([[a**2,a*b*np.cos(gam),a*c*np.cos(bet)],[a*b*np.cos(gam),b**2,b*c*np.cos(alp)],[a*c*np.cos(bet),b*c*np.cos(alp),c**2]])     
    v=d(pole1,pole2,pole3)
    
    trace_planB(pole1,pole2,pole3)
    if np.abs(alp-pi/2)<0.001 and np.abs(bet-pi/2)<0.001 and np.abs(gam-2*pi/3)<0.001:
        trace_plan(pole1,pole2,pole3)
        trace_plan(pole1,pole2,-pole3)
        trace_plan(pole2,pole1,pole3)
        trace_plan(pole2,pole1,-pole3)
        trace_plan(-pole1-pole2,pole2,pole3)
        trace_plan(-pole1-pole2,pole2,-pole3)
        trace_plan(pole1,-pole1-pole2,pole3)
        trace_plan(pole1,-pole1-pole2,-pole3)
        trace_plan(pole2,-pole1-pole2,pole3)
        trace_plan(pole2,-pole1-pole2,-pole3)
        trace_plan(-pole1-pole2,pole1,pole3)
        trace_plan(-pole1-pole2,pole1,-pole3)

    else:
        if np.abs(d(pole1,pole2,-pole3)-v)<0.001:
                trace_plan(pole1,pole2,-pole3)
        if np.abs(d(pole1,-pole2,pole3)-v)<0.001:
                trace_plan(pole1,-pole2,pole3)
        if np.abs(d(-pole1,pole2,pole3)-v)<0.001:
            trace_plan(-pole1,pole2,pole3)
        if np.abs(d(pole2,pole1,pole3)-v)<0.001:
            trace_plan(pole2,pole1,pole3)
        if np.abs(d(pole2,pole1,-pole3)-v)<0.001:
            trace_plan(pole2,pole1,-pole3)
        if np.abs(d(pole2,-pole1,pole3)-v)<0.001:
            trace_plan(pole2,-pole1,pole3)
        if np.abs(d(-pole2,pole1,pole3)-v)<0.001:
            trace_plan(-pole2,pole1,pole3)
        if np.abs(d(pole2,pole3,pole1)-v)<0.001:
            trace_plan(pole2,pole3,pole1)
        if np.abs(d(pole2,pole3,pole1)-v)<0.001:
            trace_plan(pole2,pole3,-pole1)
        if np.abs(d(pole2,-pole3,pole1)-v)<0.001:
            trace_plan(pole2,-pole3,pole1)
        if np.abs(d(-pole2,pole3,pole1)-v)<0.001:
            trace_plan(-pole2,pole3,pole1)
        if np.abs(d(pole1,pole3,pole2)-v)<0.001:
            trace_plan(pole1,pole3,pole2)
        if np.abs(d(pole1,pole3,-pole2)-v)<0.001:
            trace_plan(pole1,pole3,-pole2)
        if np.abs(d(pole1,-pole3,pole2)-v)<0.001:
            trace_plan(pole1,-pole3,pole2)
        if np.abs(d(-pole1,pole3,pole2)-v)<0.001:
            trace_plan(-pole1,pole3,pole2)
        if np.abs(d(pole3,pole1,pole2)-v)<0.001:
            trace_plan(pole3,pole1,pole2)
        if np.abs(d(pole3,pole1,-pole2)-v)<0.001:
            trace_plan(pole3,pole1,-pole2)
        if np.abs(d(pole3,-pole1,pole2)-v)<0.001:
            trace_plan(pole3,-pole1,pole2)
        if np.abs(d(-pole3,pole1,pole2)-v)<0.001:
            trace_plan(-pole3,pole1,pole2)
        if np.abs(d(pole3,pole2,pole1)-v)<0.001:
            trace_plan(pole3,pole2,pole1)
        if np.abs(d(pole3,pole2,-pole1)-v)<0.001:
            trace_plan(pole3,pole2,-pole1)
        if np.abs(d(pole3,-pole2,pole1)-v)<0.001:
            trace_plan(pole3,-pole2,pole1)
        if np.abs(d(pole3,pole2,pole1)-v)<0.001:
            trace_plan(pole3,pole2,pole1)
    trace()

def trace_plan2A(B):
    global MA,axesA,axeshA,T,VA,DA,DstarA,a
    
   
    
    for h in range(0,B.shape[0]):
        pole1=B[h,0]
        pole2=B[h,1]
        pole3=B[h,2]
        Gs=np.array([pole1,pole2,pole3],float)
        if B[h,3]==0:                    
            Gsh=np.dot(DstarA,Gs)/np.linalg.norm(np.dot(DstarA,Gs))
        else:
            Gsh=np.dot(DA,Gs)/np.linalg.norm(np.dot(DA,Gs))
        S=np.dot(MA,Gsh)
        
        if S[2]<0:
            S=-S
            Gsh=-Gsh
            pole1=-pole1
            pole2=-pole2
            pole3=-pole3
        r=np.sqrt(S[0]**2+S[1]**2+S[2]**2)
        A=np.zeros((2,100))
        Q=np.zeros((1,2))
        if S[2]==0:
             t=90
             
        else:
             t=np.arctan2(S[1],S[0])*180/pi
        w=0
        ph=np.arccos(S[2]/r)*180/pi
        for g in np.linspace(-pi,pi-0.00001,100):
            Aa=np.dot(Rot(t,0,0,1),np.dot(Rot(ph,0,1,0),np.array([np.sin(g),np.cos(g),0])))
            A[:,w]=proj(Aa[0],Aa[1],Aa[2])*600/2
            if A[0,w]<>75000:
                Q=np.vstack((Q,A[:,w]))
                w=w+1
        Q=np.delete(Q,0,0)    
        a.plot(Q[:,0]+600/2,Q[:,1]+600/2,'r')

def trace_plan2B(B):
    global MB,axesB,axeshB,T,VB,DB,DstarB,a
    
   
    
    for h in range(0,B.shape[0]):
        pole1=B[h,0]
        pole2=B[h,1]
        pole3=B[h,2]
        Gs=np.array([pole1,pole2,pole3],float)
        if B[h,3]==0:                    
            Gsh=np.dot(DstarB,Gs)/np.linalg.norm(np.dot(DstarB,Gs))
        else:
            Gsh=np.dot(DB,Gs)/np.linalg.norm(np.dot(DB,Gs))
        S=np.dot(MB,Gsh)
        
        if S[2]<0:
            S=-S
            Gsh=-Gsh
            pole1=-pole1
            pole2=-pole2
            pole3=-pole3
        r=np.sqrt(S[0]**2+S[1]**2+S[2]**2)
        A=np.zeros((2,100))
        Q=np.zeros((1,2))
        if S[2]==0:
             t=90
             
        else:
             t=np.arctan2(S[1],S[0])*180/pi
        w=0
        ph=np.arccos(S[2]/r)*180/pi
        for g in np.linspace(-pi,pi-0.00001,100):
            Aa=np.dot(Rot(t,0,0,1),np.dot(Rot(ph,0,1,0),np.array([np.sin(g),np.cos(g),0])))
            A[:,w]=proj(Aa[0],Aa[1],Aa[2])*600/2
            if A[0,w]<>75000:
                Q=np.vstack((Q,A[:,w]))
                w=w+1
        Q=np.delete(Q,0,0)    
        a.plot(Q[:,0]+600/2,Q[:,1]+600/2,'r')        

####################################################################
##### Click a pole
####################################################################    
def click_a_pole(event):
        
    global MB,DstarB,DB
    x=event.x
    y=event.y
    x=(x-411)*2/620
    y=-(y-400)*2/620
    X=2*x/(1+x**2+y**2)
    Y=2*y/(1+x**2+y**2)
    Z=(-1+x**2+y**2)/(1+x**2+y**2)
    if Z<0:
        X=-X
        Y=-Y
    A=np.dot(np.linalg.inv(MB),np.array([X,Y,Z]))
    n=0
    L=np.zeros((3,16**3))                      
                           
                        
    for i in range(-8,9,1):
        for j in range(-8,9,1):
            for k in range(-8,9,1):
                if np.linalg.norm([i,j,k])<>0:
                    if var_uvw.get()==0:
                        Q=np.dot(DstarB,np.array([i,j,k],float))/np.linalg.norm(np.dot(DstarB,np.array([i,j,k],float)))
                        if np.abs(Q[0]-A[0])<0.05 and np.abs(Q[1]-A[1])<0.05 and np.abs(Q[2]-A[2])<0.05:
                            L[:,n]=np.array([i,j,k],float)
                            n=n+1
                           
                    else:
                          
                        Q=np.dot(DB,np.array([i,j,k],float))/np.linalg.norm(np.dot(DB,np.array([i,j,k],float)))
                        if np.abs(Q[0]-A[0])<0.05 and np.abs(Q[1]-A[1])<0.05 and np.abs(Q[2]-A[2])<0.05:
                            L[:,n]=np.array([i,j,k],float)
                            n=n+1
      

    if np.linalg.norm(L[:,0])<>0:
        poleB(L[0,0],L[1,0],L[2,0])
        trace()

####################################################################
##### Fonction principale
####################################################################
def trace():
    global Ta,Tb,axesA,axeshA,MA,axesB,axeshB,MB,a,trA,trB
    a = f.add_subplot(111)
    a.figure.clear()
    a = f.add_subplot(111)
    fn = os.path.join(os.path.dirname(__file__), 'stereo.png')      
    img=Image.open(fn)
        
    Pa=np.zeros((axesA.shape[0],2))
    Pb=np.zeros((axesB.shape[0],2))
    trace_plan2A(trPA)
    trace_plan2B(trPB)
    
    for i in range(0,axesA.shape[0]):
        axeshrA=np.array([axeshA[i,0],axeshA[i,1],axeshA[i,2]])
#        print(axeshr)
        axeshrA=axeshrA/np.linalg.norm(axeshrA)
        
        Ta[i,:]=np.dot(MA,axeshrA)
        Pa[i,:]=proj(Ta[i,0],Ta[i,1],Ta[i,2])*600/2
        m=np.amax([np.abs(axesA[i,0]),np.abs(axesA[i,1]),np.abs(axesA[i,2])])
        if (np.around(axesA[i,0]/m)==axesA[i,0]/m) & (np.around(axesA[i,1]/m)==axesA[i,1]/m) & (np.around(axesA[i,2]/m)==axesA[i,2]/m):
            if var_hexaA.get()==0:
                if axeshA[i,3]==0:
                    sA='('+str(int(axesA[i,0]/m))+str(int(axesA[i,1]/m))+str(int(axesA[i,2]/m))+')'
                if axeshA[i,3]==1:
                    sA='['+str(int(axesA[i,0]/m))+str(int(axesA[i,1]/m))+str(int(axesA[i,2]/m))+']'
            if var_hexaA.get()==1:
                if axeshA[i,3]==0:
                    sA='('+str(int(axesA[i,0]/m))+str(int(axesA[i,1]/m))+str(-int(axesA[i,1]/m)-int(axesA[i,0]/m))+str(int(axesA[i,2]/m))+')'
                if axeshA[i,3]==1:
                    sA='['+str(int(2*(axesA[i,0]/m)-axesA[i,1]/m))+str(int(2*(axesA[i,1]/m)-axesA[i,0]/m))+str(-int(axesA[i,1]/m)-int(axesA[i,0]/m))+str(int(3*axesA[i,2]/m))+']'
                
        else:
            if var_hexaA.get()==0:
                if axeshA[i,3]==0:
                    sA='('+str(int(axesA[i,0]))+str(int(axesA[i,1]))+str(int(axesA[i,2]))+')'
                if axeshA[i,3]==1:
                    sA='['+str(int(axesA[i,0]))+str(int(axesA[i,1]))+str(int(axesA[i,2]))+']'
            if var_hexaA.get()==1:
                if axeshA[i,3]==0:
                    sA='('+str(int(axesA[i,0]))+str(int(axesA[i,1]))+str(-int(axesA[i,1])-int(axesA[i,0]))+str(int(axesA[i,2]))+')'
                if axeshA[i,3]==1:
                    sA='['+str(int(2*(axesA[i,0])-axesA[i,1]))+str(int(2*(axesA[i,1])-axesA[i,0]))+str(-int(axesA[i,1])-int(axesA[i,0]))+str(int(3*axesA[i,2]))+']'
                
       
       
        a.annotate(sA,(Pa[i,0]+610/2,Pa[i,1]+610/2))
    for i in range(0,axesB.shape[0]):
        axeshrB=np.array([axeshB[i,0],axeshB[i,1],axeshB[i,2]])
#        print(axeshr)
        axeshrB=axeshrB/np.linalg.norm(axeshrB)
        
        Tb[i,:]=np.dot(MB,axeshrB)
        Pb[i,:]=proj(Tb[i,0],Tb[i,1],Tb[i,2])*600/2
        m=np.amax([np.abs(axesB[i,0]),np.abs(axesB[i,1]),np.abs(axesB[i,2])])
        if (np.around(axesB[i,0]/m)==axesB[i,0]/m) & (np.around(axesB[i,1]/m)==axesB[i,1]/m) & (np.around(axesB[i,2]/m)==axesB[i,2]/m):
            if var_hexaB.get()==0:
                if axeshB[i,3]==0:
                    sB='('+str(int(axesB[i,0]/m))+str(int(axesB[i,1]/m))+str(int(axesB[i,2]/m))+')'
                if axeshB[i,3]==1:
                    sB='['+str(int(axesB[i,0]/m))+str(int(axesB[i,1]/m))+str(int(axesB[i,2]/m))+']'
            if var_hexaB.get()==1:
                if axeshB[i,3]==0:
                    sB='('+str(int(axesB[i,0]/m))+str(int(axesB[i,1]/m))+str(-int(axesB[i,1]/m)-int(axesB[i,0]/m))+str(int(axesB[i,2]/m))+')'
                if axeshB[i,3]==1:
                    sB='['+str(int(2*(axesB[i,0]/m)-axesB[i,1]/m))+str(int(2*(axesB[i,1]/m)-axesB[i,0]/m))+str(-int(axesB[i,1]/m)-int(axesB[i,0]/m))+str(int(3*axesB[i,2]/m))+']'
                
        else:
            if var_hexaB.get()==0:
                if axeshB[i,3]==0:
                    sB='('+str(int(axesB[i,0]))+str(int(axesB[i,1]))+str(int(axesB[i,2]))+')'
                if axeshB[i,3]==1:
                    sB='['+str(int(axesB[i,0]))+str(int(axesB[i,1]))+str(int(axesB[i,2]))+']'
            if var_hexaB.get()==1:
                if axeshB[i,3]==0:
                    sB='('+str(int(axesB[i,0]))+str(int(axesB[i,1]))+str(-int(axesB[i,1])-int(axesB[i,0]))+str(int(axesB[i,2]))+')'
                if axeshB[i,3]==1:
                    sB='['+str(int(2*(axesB[i,0])-axesB[i,1]))+str(int(2*(axesB[i,1])-axesB[i,0]))+str(-int(axesB[i,1])-int(axesB[i,0]))+str(int(3*axesB[i,2]))+']'
                
       
        a.annotate(sB,(Pb[i,0]+610/2,Pb[i,1]+610/2-20))
                
        
    a.plot(Pa[:,0]+600/2,Pa[:,1]+600/2,'bo', ms=10)
    a.plot(Pb[:,0]+600/2,Pb[:,1]+600/2,'bo', mfc='none', mec='b', markeredgewidth='1.5',ms=10)
    
    a.axis([0,600,0,600])
    a.imshow(img)
    a.axis('off')   
    a.figure.canvas.draw()
    

    
def princ():
    global Ta,Tb,MA,MB,trPA,trPB
    trPA=np.zeros((1,4))
    trPB=np.zeros((1,4))
    a = f.add_subplot(111)
    a.figure.clear()
    a = f.add_subplot(111)
    phi1a=eval(phi1A_entry.get())
    phia=eval(phiA_entry.get())
    phi2a=eval(phi2A_entry.get())
    phi1b=eval(phi1B_entry.get())
    phib=eval(phiB_entry.get())
    phi2b=eval(phi2B_entry.get())
    fn = os.path.join(os.path.dirname(__file__), 'stereo.png')      
    img=Image.open(fn)
    cristA()    
    cristB()
    Pa=np.zeros((axesA.shape[0],2))
    Ta=np.zeros((axesA.shape))
    Pb=np.zeros((axesB.shape[0],2))
    Tb=np.zeros((axesB.shape))
    
    for i in range(0,axesA.shape[0]):
        axeshrA=np.array([axeshA[i,0],axeshA[i,1],axeshA[i,2]])
#        print(axeshr)
        axeshrA=axeshrA/np.linalg.norm(axeshrA)
        
        Ta[i,:]=np.dot(rotation(phi1a,phia,phi2a),axeshrA)
        Pa[i,:]=proj(Ta[i,0],Ta[i,1],Ta[i,2])*600/2
        m=np.amax([np.abs(axesA[i,0]),np.abs(axesA[i,1]),np.abs(axesA[i,2])])
        if (np.around(axesA[i,0]/m)==axesA[i,0]/m) & (np.around(axesA[i,1]/m)==axesA[i,1]/m) & (np.around(axesA[i,2]/m)==axesA[i,2]/m):
            if var_hexaA.get()==0:
                if axeshA[i,3]==0:
                    sA='('+str(int(axesA[i,0]/m))+str(int(axesA[i,1]/m))+str(int(axesA[i,2]/m))+')'
                if axeshA[i,3]==1:
                    sA='['+str(int(axesA[i,0]/m))+str(int(axesA[i,1]/m))+str(int(axesA[i,2]/m))+']'
            if var_hexaA.get()==1:
                if axeshA[i,3]==0:
                    sA='('+str(int(axesA[i,0]/m))+str(int(axesA[i,1]/m))+str(-int(axesA[i,1]/m)-int(axesA[i,0]/m))+str(int(axesA[i,2]/m))+')'
                if axeshA[i,3]==1:
                    sA='['+str(int(2*(axesA[i,0]/m)-axesA[i,1]/m))+str(int(2*(axesA[i,1]/m)-axesA[i,0]/m))+str(-int(axesA[i,1]/m)-int(axesA[i,0]/m))+str(int(3*axesA[i,2]/m))+']'
                
        else:
            if var_hexaA.get()==0:
                if axeshA[i,3]==0:
                    sA='('+str(int(axesA[i,0]))+str(int(axesA[i,1]))+str(int(axesA[i,2]))+')'
                if axeshA[i,3]==1:
                    sA='['+str(int(axesA[i,0]))+str(int(axesA[i,1]))+str(int(axesA[i,2]))+']'
            if var_hexaA.get()==1:
                if axeshA[i,3]==0:
                    sA='('+str(int(axesA[i,0]))+str(int(axesA[i,1]))+str(-int(axesA[i,1])-int(axesA[i,0]))+str(int(axesA[i,2]))+')'
                if axeshA[i,3]==1:
                    sA='['+str(int(2*(axesA[i,0])-axesA[i,1]))+str(int(2*(axesA[i,1])-axesA[i,0]))+str(-int(axesA[i,1])-int(axesA[i,0]))+str(int(3*axesA[i,2]))+']'
                
       
        a.annotate(sA,(Pa[i,0]+610/2,Pa[i,1]+610/2))
    for i in range(0,axesB.shape[0]):
        axeshrB=np.array([axeshB[i,0],axeshB[i,1],axeshB[i,2]])
#        print(axeshr)
        axeshrB=axeshrB/np.linalg.norm(axeshrB)
        
        Tb[i,:]=np.dot(rotation(phi1b,phib,phi2b),axeshrB)
        Pb[i,:]=proj(Tb[i,0],Tb[i,1],Tb[i,2])*600/2
        m=np.amax([np.abs(axesB[i,0]),np.abs(axesB[i,1]),np.abs(axesB[i,2])])
        if (np.around(axesB[i,0]/m)==axesB[i,0]/m) & (np.around(axesB[i,1]/m)==axesB[i,1]/m) & (np.around(axesB[i,2]/m)==axesB[i,2]/m):
            if var_hexaB.get()==0:
                if axeshB[i,3]==0:
                    sB='('+str(int(axesB[i,0]/m))+str(int(axesB[i,1]/m))+str(int(axesB[i,2]/m))+')'
                if axeshB[i,3]==1:
                    sB='['+str(int(axesB[i,0]/m))+str(int(axesB[i,1]/m))+str(int(axesB[i,2]/m))+']'
            if var_hexaB.get()==1:
                if axeshB[i,3]==0:
                    sB='('+str(int(axesB[i,0]/m))+str(int(axesB[i,1]/m))+str(-int(axesB[i,1]/m)-int(axesB[i,0]/m))+str(int(axesB[i,2]/m))+')'
                if axeshB[i,3]==1:
                    sB='['+str(int(2*(axesB[i,0]/m)-axesB[i,1]/m))+str(int(2*(axesB[i,1]/m)-axesB[i,0]/m))+str(-int(axesB[i,1]/m)-int(axesB[i,0]/m))+str(int(3*axesB[i,2]/m))+']'
                
        else:
            if var_hexaB.get()==0:
                if axeshB[i,3]==0:
                    sB='('+str(int(axesB[i,0]))+str(int(axesB[i,1]))+str(int(axesB[i,2]))+')'
                if axeshB[i,3]==1:
                    sB='['+str(int(axesB[i,0]))+str(int(axesB[i,1]))+str(int(axesB[i,2]))+']'
            if var_hexaB.get()==1:
                if axeshB[i,3]==0:
                    sB='('+str(int(axesB[i,0]))+str(int(axesB[i,1]))+str(-int(axesB[i,1])-int(axesB[i,0]))+str(int(axesB[i,2]))+')'
                if axeshB[i,3]==1:
                    sB='['+str(int(2*(axesB[i,0])-axesB[i,1]))+str(int(2*(axesB[i,1])-axesB[i,0]))+str(-int(axesB[i,1])-int(axesB[i,0]))+str(int(3*axesB[i,2]))+']'
                
       
        a.annotate(sB,(Pb[i,0]+610/2,Pb[i,1]+610/2-20))
                
               
    a.plot(Pa[:,0]+600/2,Pa[:,1]+600/2,'bo', ms=10)
    a.plot(Pb[:,0]+600/2,Pb[:,1]+600/2,'bo', mfc='none', mec='b', markeredgewidth='1.5',ms=10)
    a.axis([0,600,0,600])
    
    a.imshow(img)
    a.axis('off')   
    a.figure.canvas.draw() 
    
        
    MA=rotation(phi1a,phia,phi2a)
    MB=rotation(phi1b,phib,phi2b)
   
    return Ta,MA,MB,Tb
    
                    
######################################################################
# GUI
######################################################################

#def file_save():
#    global D1,D0
#    fout = asksaveasfile(mode='w', defaultextension=".txt")
#    
#    for i in range(D1.shape[0]):
#        text2save = str(int(D0[i,4]))+'\t'+'['+str(int(D1[i,0]))+','+str(int(D1[i,1]))+','+str(int(D1[i,2]))+']'+'\t '+str(np.around(D0[i,3],decimals=2))
#        fout.write("%s\n" % text2save)
#        
#    fout.close()    

def image_save():
    
    s = asksaveasfile(mode='w', defaultextension=".jpg")    
    if s:    
        f.savefig(s.name)
    #s.close()
    
####################################################
#fonction d'initialisation
##################################################
def init():
    global var_uvw,dmipA,dmipB,dA_label_var,dB_label_var,var_hexaA,var_hexaB,trPA, trPB
    fn = os.path.join(os.path.dirname(__file__), 'stereo.png')      
    img=Image.open(fn)
    a = f.add_subplot(111)
    a.axis('off')
    a.imshow(img)
    a.figure.canvas.draw()
    
    var_uvw=IntVar()
    trPA=np.zeros((1,3))
    trPB=np.zeros((1,3))
    dA_label_var=StringVar()
    dB_label_var=StringVar()
    
    dA_label_var.set(0)
    dB_label_var.set(0)
    var_hexaA=IntVar()
    var_hexaB=IntVar()
    dmipA=0
    dmipB=0
    return var_uvw
   

##############################################################
# fonction pour quitter
#######################################################
def _quit():
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate
#############################################################


root = Tk()
root.wm_title("Stereo-inter")
root.geometry('1220x798+10+40')
root.configure(bg = '#BDBDBD')
#root.resizable(0,0)
#s=ttk.Style()
#s.theme_use('clam')
style = ttk.Style()
theme = style.theme_use()
default = style.lookup(theme, 'background')



################################################
# Creation d'une zone pour tracer des graphiques
################################################
f = Figure(facecolor='white',figsize=[2,2],dpi=100)

canvas = FigureCanvasTkAgg(f, master=root)

canvas.get_tk_widget().place(x=0,y=0,height=800,width=800)
canvas._tkcanvas.bind('<Button-3>', click_a_pole)
canvas.show()
#toolbar = NavigationToolbar2TkAgg( canvas, root )
#toolbar.zoom('off')
#toolbar.update()


###################################################

init()
#import _imaging
#print _imaging.__file__

##############################################
# Boutons
##############################################


label_euler = Label (master=root)
label_euler.place(relx=0.77,rely=0.45,height=19,width=75)
label_euler.configure(activebackground="#cccccc")
label_euler.configure(activeforeground="black")
label_euler.configure(cursor="fleur")
label_euler.configure(foreground="black")
label_euler.configure(highlightcolor="black")
label_euler.configure(text='''Euler angles''')



phi1A_entry = Entry (master=root)
phi1A_entry.place(relx=0.7,rely=0.5,relheight=0.03,relwidth=0.07)
phi1A_entry.configure(background="white")
phi1A_entry.configure(foreground="black")
phi1A_entry.configure(highlightbackground="#e0e0dfdfe3e3")
phi1A_entry.configure(highlightcolor="#000000")
phi1A_entry.configure(insertbackground="#000000")
phi1A_entry.configure(selectbackground="#c4c4c4")
phi1A_entry.configure(selectforeground="black")

phiA_entry = Entry (master=root)
phiA_entry.place(relx=0.7,rely=0.55,relheight=0.03,relwidth=0.07)
phiA_entry.configure(background="white")
phiA_entry.configure(foreground="black")
phiA_entry.configure(highlightcolor="black")
phiA_entry.configure(insertbackground="black")
phiA_entry.configure(selectbackground="#c4c4c4")
phiA_entry.configure(selectforeground="black")

phi2A_entry = Entry (master=root)
phi2A_entry.place(relx=0.7,rely=0.6,relheight=0.03,relwidth=0.07)
phi2A_entry.configure(background="white")
phi2A_entry.configure(foreground="black")
phi2A_entry.configure(highlightcolor="black")
phi2A_entry.configure(insertbackground="black")
phi2A_entry.configure(selectbackground="#c4c4c4")
phi2A_entry.configure(selectforeground="black")


Phi1A_label = Label (master=root)
Phi1A_label.place(relx=0.67,rely=0.5,height=19,width=33)
Phi1A_label.configure(activebackground="#cccccc")
Phi1A_label.configure(activeforeground="black")
Phi1A_label.configure(foreground="black")
Phi1A_label.configure(highlightcolor="black")
Phi1A_label.configure(text='''Phi1A''')

PhiA_label = Label (master=root)
PhiA_label.place(relx=0.67,rely=0.55,height=19,width=27)
PhiA_label.configure(activebackground="#cccccc")
PhiA_label.configure(activeforeground="black")
PhiA_label.configure(foreground="black")
PhiA_label.configure(highlightcolor="black")
PhiA_label.configure(text='''PhiA''')

Phi2A_label = Label (master=root)
Phi2A_label.place(relx=0.67,rely=0.6,height=19,width=33)
Phi2A_label.configure(activebackground="#cccccc")
Phi2A_label.configure(activeforeground="black")
Phi2A_label.configure(foreground="black")
Phi2A_label.configure(highlightcolor="black")
Phi2A_label.configure(text='''Phi2A''')
phi1B_entry = Entry (master=root)
phi1B_entry.place(relx=0.84,rely=0.5,relheight=0.03,relwidth=0.07)
phi1B_entry.configure(background="white")

phi1B_entry.configure(foreground="black")
phi1B_entry.configure(highlightbackground="#e0e0dfdfe3e3")
phi1B_entry.configure(highlightcolor="#000000")
phi1B_entry.configure(insertbackground="#000000")
phi1B_entry.configure(selectbackground="#c4c4c4")
phi1B_entry.configure(selectforeground="black")

Phi1B = Label (master=root)
Phi1B.place(relx=0.8,rely=0.5,height=29,width=33)
Phi1B.configure(activebackground="#cccccc")
Phi1B.configure(activeforeground="black")
Phi1B.configure(foreground="black")
Phi1B.configure(highlightcolor="black")
Phi1B.configure(text='''Phi1B''')

PhiB_label1 = Label (master=root)
PhiB_label1.place(relx=0.81,rely=0.55,height=19,width=26)
PhiB_label1.configure(activebackground="#cccccc")
PhiB_label1.configure(activeforeground="black")
PhiB_label1.configure(foreground="black")
PhiB_label1.configure(highlightcolor="black")
PhiB_label1.configure(text='''PhiB''')

Phi2B_label2 = Label (master=root)
Phi2B_label2.place(relx=0.8,rely=0.6,height=19,width=32)
Phi2B_label2.configure(activebackground="#cccccc")
Phi2B_label2.configure(activeforeground="black")
Phi2B_label2.configure(foreground="black")
Phi2B_label2.configure(highlightcolor="black")
Phi2B_label2.configure(text='''Phi2B''')

phiB_entry = Entry (master=root)
phiB_entry.place(relx=0.84,rely=0.55,relheight=0.03,relwidth=0.07)
phiB_entry.configure(background="white")

phiB_entry.configure(foreground="black")
phiB_entry.configure(highlightbackground="#e0e0dfdfe3e3")
phiB_entry.configure(highlightcolor="#000000")
phiB_entry.configure(insertbackground="#000000")
phiB_entry.configure(selectbackground="#c4c4c4")
phiB_entry.configure(selectforeground="black")

phi2B_entry = Entry (master=root)
phi2B_entry.place(relx=0.84,rely=0.6,relheight=0.03,relwidth=0.07)
phi2B_entry.configure(background="white")

phi2B_entry.configure(foreground="black")
phi2B_entry.configure(highlightbackground="#e0e0dfdfe3e3")
phi2B_entry.configure(highlightcolor="#000000")
phi2B_entry.configure(insertbackground="#000000")
phi2B_entry.configure(selectbackground="#c4c4c4")
phi2B_entry.configure(selectforeground="black")

button_trace = Button (master=root)
button_trace.place(relx=0.7,rely=0.66,height=21,width=49)
button_trace.configure(activebackground="#f9f9f9")
button_trace.configure(activeforeground="black")
button_trace.configure(background="#ff0000")
button_trace.configure(command=princ)
button_trace.configure(foreground="black")
button_trace.configure(highlightcolor="black")
button_trace.configure(pady="0")
button_trace.configure(text='''Plot''')
Cristal_label = Label (master=root)
Cristal_label.place(relx=0.66,rely=0.03,height=19,width=142)
Cristal_label.configure(text='''Crystal parameters''')

a_cristal_label = Label (master=root)
a_cristal_label.place(relx=0.68,rely=0.06,height=19,width=12)
a_cristal_label.configure(text='''a''')

b_cristal_label = Label (master=root)
b_cristal_label.place(relx=0.68,rely=0.1,height=19,width=12)
b_cristal_label.configure(activebackground="#f9f9f9")
b_cristal_label.configure(activeforeground="black")
b_cristal_label.configure(foreground="black")
b_cristal_label.configure(highlightcolor="black")
b_cristal_label.configure(text='''b''')

c_cristal_label = Label (master=root)
c_cristal_label.place(relx=0.68,rely=0.14,height=19,width=11)
c_cristal_label.configure(activebackground="#f9f9f9")
c_cristal_label.configure(activeforeground="black")
c_cristal_label.configure(foreground="black")
c_cristal_label.configure(highlightcolor="black")
c_cristal_label.configure(text='''c''')

alp_cristal_label = Label (master=root)
alp_cristal_label.place(relx=0.67,rely=0.19,height=19,width=35)
alp_cristal_label.configure(activebackground="#f9f9f9")
alp_cristal_label.configure(activeforeground="black")
alp_cristal_label.configure(foreground="black")
alp_cristal_label.configure(highlightcolor="black")
alp_cristal_label.configure(text='''alpha''')

bet_cristal_label = Label (master=root)
bet_cristal_label.place(relx=0.67,rely=0.23,height=19,width=32)
bet_cristal_label.configure(activebackground="#f9f9f9")
bet_cristal_label.configure(activeforeground="black")
bet_cristal_label.configure(foreground="black")
bet_cristal_label.configure(highlightcolor="black")
bet_cristal_label.configure(text='''beta''')

gam_cristal_label = Label (master=root)
gam_cristal_label.place(relx=0.66,rely=0.26,height=19,width=45)
gam_cristal_label.configure(activebackground="#f9f9f9")
gam_cristal_label.configure(activeforeground="black")
gam_cristal_label.configure(foreground="black")
gam_cristal_label.configure(highlightcolor="black")
gam_cristal_label.configure(text='''gamma''')

aA_entry = Entry (master=root)
aA_entry.place(relx=0.7,rely=0.06,relheight=0.03,relwidth=0.03)

aB_entry = Entry (master=root)
aB_entry.place(relx=0.75,rely=0.06,relheight=0.03,relwidth=0.03)

bA_entry = Entry (master=root)
bA_entry.place(relx=0.7,rely=0.1,relheight=0.03,relwidth=0.03)

bB_entry = Entry (master=root)
bB_entry.place(relx=0.75,rely=0.1,relheight=0.03,relwidth=0.03)


cA_entry = Entry (master=root)
cA_entry.place(relx=0.7,rely=0.14,relheight=0.03,relwidth=0.03)

cB_entry = Entry (master=root)
cB_entry.place(relx=0.75,rely=0.14,relheight=0.03,relwidth=0.03)


alpA_entry = Entry (master=root)
alpA_entry.place(relx=0.7,rely=0.18,relheight=0.03,relwidth=0.03)

alpB_entry = Entry (master=root)
alpB_entry.place(relx=0.75,rely=0.18,relheight=0.03,relwidth=0.03)

betA_entry = Entry (master=root)
betA_entry.place(relx=0.7,rely=0.22,relheight=0.03,relwidth=0.03)

betB_entry = Entry (master=root)
betB_entry.place(relx=0.75,rely=0.22,relheight=0.03,relwidth=0.03)

gamA_entry = Entry (master=root)
gamA_entry.place(relx=0.7,rely=0.26,relheight=0.03,relwidth=0.03)

gamB_entry = Entry (master=root)
gamB_entry.place(relx=0.75,rely=0.26,relheight=0.03,relwidth=0.03)

uvw_button = Checkbutton (master=root)
uvw_button.place(relx=0.75,rely=0.66,relheight=0.03,relwidth=0.04)
uvw_button.configure(text='''uvw''')
uvw_button.configure(variable=var_uvw)

e_label = Label (master=root)
e_label.place(relx=0.66,rely=0.31,height=19,width=56)
e_label.configure(text='''Max indices''')

eA_entry = Entry (master=root)
eA_entry.place(relx=0.72,rely=0.31,relheight=0.03,relwidth=0.02)


eB_entry = Entry (master=root)
eB_entry.place(relx=0.77,rely=0.31,relheight=0.03,relwidth=0.02)

e2_label = Label (master=root)
e2_label.place(relx=0.68,rely=0.36,height=19,width=12)
e2_label.configure(text='''d''')

dmA_button = Button (master=root)
dmA_button.place(relx=0.7,rely=0.36,height=21,width=13)
dmA_button.configure(command=dmA)
dmA_button.configure(pady="0")
dmA_button.configure(text='''-''')

dmB_button = Button (master=root)
dmB_button.place(relx=0.8,rely=0.36,height=21,width=13)
dmB_button.configure(command=dmB)
dmB_button.configure(pady="0")
dmB_button.configure(text='''-''')

dA_entry = Entry (master=root)
dA_entry.place(relx=0.72,rely=0.36,relheight=0.02,relwidth=0.04)

dB_entry = Entry (master=root)
dB_entry.place(relx=0.82,rely=0.36,relheight=0.02,relwidth=0.04)

dpA_button = Button (master=root)
dpA_button.place(relx=0.76,rely=0.36,height=21,width=17)
dpA_button.configure(command=dpA)
dpA_button.configure(pady="0")
dpA_button.configure(text='''+''')

dpB_button = Button (master=root)
dpB_button.place(relx=0.86,rely=0.36,height=21,width=17)
dpB_button.configure(command=dpB)
dpB_button.configure(pady="0")
dpB_button.configure(text='''+''')

dA_label = Label (master=root)
dA_label.place(relx=0.73,rely=0.39,height=19,width=16)
dA_label.configure(textvariable=dA_label_var)

dB_label = Label (master=root)
dB_label.place(relx=0.83,rely=0.39,height=19,width=16)
dB_label.configure(textvariable=dB_label_var)


label_addpoleA = Label (master=root)
label_addpoleA.place(relx=0.81,rely=0.03,height=19,width=90)
label_addpoleA.configure(activebackground="#cccccc")
label_addpoleA.configure(activeforeground="black")
label_addpoleA.configure(foreground="black")
label_addpoleA.configure(highlightcolor="black")
label_addpoleA.configure(text='''Add pole A''')

pole1A_entry = Entry (master=root)
pole1A_entry.place(relx=0.81,rely=0.06,relheight=0.02
,relwidth=0.04)
pole1A_entry.configure(background="white")
pole1A_entry.configure(foreground="black")
pole1A_entry.configure(highlightcolor="black")
pole1A_entry.configure(insertbackground="black")
pole1A_entry.configure(selectbackground="#c4c4c4")
pole1A_entry.configure(selectforeground="black")

pole2A_entry = Entry (master=root)
pole2A_entry.place(relx=0.87,rely=0.06,relheight=0.02
,relwidth=0.04)
pole2A_entry.configure(background="white")
pole2A_entry.configure(foreground="black")
pole2A_entry.configure(highlightcolor="black")
pole2A_entry.configure(insertbackground="black")
pole2A_entry.configure(selectbackground="#c4c4c4")
pole2A_entry.configure(selectforeground="black")

pole3A_entry = Entry (master=root)
pole3A_entry.place(relx=0.93,rely=0.06,relheight=0.02
,relwidth=0.04)
pole3A_entry.configure(background="white")
pole3A_entry.configure(foreground="black")
pole3A_entry.configure(highlightcolor="black")
pole3A_entry.configure(insertbackground="black")
pole3A_entry.configure(selectbackground="#c4c4c4")
pole3A_entry.configure(selectforeground="black")

addpoleA_button = Button (master=root)
addpoleA_button.place(relx=0.81,rely=0.11,height=31,width=57)
addpoleA_button.configure(activebackground="#f9f9f9")
addpoleA_button.configure(activeforeground="black")
addpoleA_button.configure(command=addpoleA)
addpoleA_button.configure(foreground="black")
addpoleA_button.configure(highlightcolor="black")
addpoleA_button.configure(pady="0")
addpoleA_button.configure(text='''Add''')

symA_button = Button (master=root)
symA_button.place(relx=0.87,rely=0.11,height=31,width=71)
symA_button.configure(command=addpoleA_sym)
symA_button.configure(pady="0")
symA_button.configure(text='''Symmetry''')

trace_planA_button = Button (master=root)
trace_planA_button.place(relx=0.93,rely=0.11,height=31,width=61)
trace_planA_button.configure(command=trace_addplanA)
trace_planA_button.configure(pady="0")
trace_planA_button.configure(text='''Plane''')


label_addpoleB = Label (master=root)
label_addpoleB.place(relx=0.81,rely=0.2,height=19,width=90)
label_addpoleB.configure(activebackground="#cccccc")
label_addpoleB.configure(activeforeground="black")
label_addpoleB.configure(foreground="black")
label_addpoleB.configure(highlightcolor="black")
label_addpoleB.configure(text='''Add pole B''')

pole1B_entry = Entry (master=root)
pole1B_entry.place(relx=0.81,rely=0.24,relheight=0.02
,relwidth=0.04)
pole1B_entry.configure(background="white")
pole1B_entry.configure(foreground="black")
pole1B_entry.configure(highlightcolor="black")
pole1B_entry.configure(insertbackground="black")
pole1B_entry.configure(selectbackground="#c4c4c4")
pole1B_entry.configure(selectforeground="black")

pole2B_entry = Entry (master=root)
pole2B_entry.place(relx=0.87,rely=0.24,relheight=0.02
,relwidth=0.04)
pole2B_entry.configure(background="white")
pole2B_entry.configure(foreground="black")
pole2B_entry.configure(highlightcolor="black")
pole2B_entry.configure(insertbackground="black")
pole2B_entry.configure(selectbackground="#c4c4c4")
pole2B_entry.configure(selectforeground="black")

pole3B_entry = Entry (master=root)
pole3B_entry.place(relx=0.93,rely=0.24,relheight=0.02
,relwidth=0.04)
pole3B_entry.configure(background="white")
pole3B_entry.configure(foreground="black")
pole3B_entry.configure(highlightcolor="black")
pole3B_entry.configure(insertbackground="black")
pole3B_entry.configure(selectbackground="#c4c4c4")
pole3B_entry.configure(selectforeground="black")

addpoleB_button = Button (master=root)
addpoleB_button.place(relx=0.81,rely=0.28,height=31,width=55)
addpoleB_button.configure(activebackground="#f9f9f9")
addpoleB_button.configure(activeforeground="black")
addpoleB_button.configure(command=addpoleB)
addpoleB_button.configure(foreground="black")
addpoleB_button.configure(highlightcolor="black")
addpoleB_button.configure(pady="0")
addpoleB_button.configure(text='''Add''')

symB_button = Button (master=root)
symB_button.place(relx=0.87,rely=0.28,height=31,width=67)
symB_button.configure(command=addpoleB_sym)
symB_button.configure(pady="0")
symB_button.configure(text='''Symmetry''')

trace_planB_button = Button (master=root)
trace_planB_button.place(relx=0.93,rely=0.28,height=31,width=59)
trace_planB_button.configure(command=trace_addplanB)
trace_planB_button.configure(pady="0")
trace_planB_button.configure(text='''Plane''')

hexaA_button = Checkbutton (master=root)
hexaA_button.place(relx=0.73,rely=0.71,relheight=0.03,relwidth=0.06)
hexaA_button.configure(text='''Hexa A''')
hexaA_button.configure(variable=var_hexaA)

hexaB_button = Checkbutton (master=root)
hexaB_button.place(relx=0.79,rely=0.71,relheight=0.03,relwidth=0.06)
hexaB_button.configure(text='''Hexa B''')
hexaB_button.configure(variable=var_hexaB)

menu = Menu(master=root)
filemenu = Menu(menu, tearoff=0)
menu.add_cascade(label="Save", menu=filemenu)

root.config(menu=menu)

filemenu.add_command(label="Save figure", command=image_save) 
######################################################################################################
######## importer des structures cristallines depuis un fichier Nom,a,b,c,alpha,beta,gamma,space group
######################################################################################################

def structureA(i0):
    global x0
    
    aA_entry.delete(0,END)
    aA_entry.insert(1,eval(x0[i0][1]))
    bA_entry.delete(0,END)    
    bA_entry.insert(1,eval(x0[i0][2]))
    cA_entry.delete(0,END)    
    cA_entry.insert(1,eval(x0[i0][3]))
    alpA_entry.delete(0,END)    
    alpA_entry.insert(1,eval(x0[i0][4]))
    betA_entry.delete(0,END)    
    betA_entry.insert(1,eval(x0[i0][5]))
    gamA_entry.delete(0,END)    
    gamA_entry.insert(1,eval(x0[i0][6]))

def structureB(i0):
    global x0
    
    aB_entry.delete(0,END)
    aB_entry.insert(1,eval(x0[i0][1]))
    bB_entry.delete(0,END)    
    bB_entry.insert(1,eval(x0[i0][2]))
    cB_entry.delete(0,END)    
    cB_entry.insert(1,eval(x0[i0][3]))
    alpB_entry.delete(0,END)    
    alpB_entry.insert(1,eval(x0[i0][4]))
    betB_entry.delete(0,END)    
    betB_entry.insert(1,eval(x0[i0][5]))
    gamB_entry.delete(0,END)    
    gamB_entry.insert(1,eval(x0[i0][6]))
    
def createstructureA(i):
    return lambda:structureA(i)

def createstructureB(i):
    return lambda:structureB(i)  
  
cristalmenuA=Menu(menu,tearoff=0)  
cristalmenuB=Menu(menu,tearoff=0)
menu.add_cascade(label="Structures A", menu=cristalmenuA)
menu.add_cascade(label="Structures B", menu=cristalmenuB)
file_struct=open(os.path.join(os.path.dirname(__file__), 'structure.txt') ,"r")

x0=[]

i=0
for line in file_struct:
    x0.append(map(str, line.split()))
    cristalmenuA.add_command(label=x0[i][0], command=createstructureA(i))
    cristalmenuB.add_command(label=x0[i][0], command=createstructureB(i))
    i=i+1

file_struct.close()

#######################################################################################################

phi1A_entry.insert(0,0)
phiA_entry.insert(0,0)
phi2A_entry.insert(0,0)
phi1B_entry.insert(0,0)
phiB_entry.insert(0,0)
phi2B_entry.insert(0,0)
eA_entry.insert(1,1)
dA_entry.insert(1,1)
eB_entry.insert(1,1)
dB_entry.insert(1,1)

mainloop()     
              
              
              
