#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class for creating horseshoe vortex panels to be used in a vortex lattice method

Author: David Darmofal, MIT
Date: April 2, 2021

"""
import numpy as np

class HorseshoePanel(object):
    def __init__(self, XLEa, XLEb, ca, cb, aga, agb, al0a, al0b):
        self.XLEa = np.array(XLEa, dtype=np.float64)
        self.XLEb = np.array(XLEb, dtype=np.float64)
        self.ca   = ca
        self.cb   = cb
        self.aga  = aga
        self.agb  = agb
        self.al0a = al0a
        self.al0b = al0b
        
    def get_XLEa(self):
        return self.XLEa.copy()
    
    def get_XLEb(self):
        return self.XLEb.copy()
    
    def get_ca(self):
        return self.ca
    
    def get_cb(self):
        return self.cb   

    def get_aga(self):
        return self.aga
    
    def get_agb(self):
        return self.agb      

    def get_al0a(self):
        return self.al0a
    
    def get_al0b(self):
        return self.al0b  
    
    def get_Xcontrol(self):
        Xcontrol = 0.5*(self.XLEa+self.XLEb)
        Xcontrol[0] += 0.375*(self.ca+self.cb)
        return Xcontrol
    
    def get_ncontrol(self):
        ag  = 0.5*(self.aga + self.agb)
        al0 = 0.5*(self.al0a + self.al0b)
        ac = (ag-al0)*np.pi/180.
        nax = ac # small angle approximation
        naz = 1.0 # small angle approximation
        dy = self.XLEb[1]-self.XLEa[1]
        dz = self.XLEb[2]-self.XLEa[2]
        lyz = (dy**2 + dz**2)**(0.5)
        cth = dy/lyz
        sth = dz/lyz
        ncx = nax
        ncy = -sth*naz
        ncz =  cth*naz
        return np.array([ncx, ncy, ncz])
    
    def get_Xbound(self):
        Xbound = 0.5*(self.XLEa+self.XLEb)
        Xbound[0] += 0.125*(self.ca+self.cb)
        return Xbound
    
    def get_Xbounda(self):
        Xbounda = self.XLEa.copy()
        Xbounda[0] += 0.25*self.ca
        return Xbounda
    
    def get_Xboundb(self):
        Xboundb = self.XLEb.copy()
        Xboundb[0] += 0.25*self.cb
        return Xboundb
    
    def get_unitV(self, Xc):
        Xa = self.get_Xbounda()
        Xb = self.get_Xboundb()
        a = Xc - Xa
        b = Xc - Xb
        x = np.array([1.,0,0])
        
        amag = np.linalg.norm(a)
        bmag = np.linalg.norm(b)
        
        Vax =  np.cross(a,x)/(amag-a[0])/amag
        Vbx = -np.cross(b,x)/(bmag-b[0])/bmag
        Vab =  np.cross(a,b)/(amag*bmag+np.dot(a,b))*(1/amag + 1/bmag)
        
        unitV = (Vax + Vbx + Vab)/(4*np.pi)
    
        return unitV
    
    def get_unitVind(self, Xc):
        Xa = self.get_Xbounda()
        Xb = self.get_Xboundb()
        a = Xc - Xa
        b = Xc - Xb
        x = np.array([1.,0,0])
        
        amag = np.linalg.norm(a)
        bmag = np.linalg.norm(b)
        
        Vax =  np.cross(a,x)/(amag-a[0])/amag
        Vbx = -np.cross(b,x)/(bmag-b[0])/bmag
        
        unitVind = (Vax + Vbx)/(4*np.pi)
    
        return unitVind
        
if __name__ == "__main__":
    H = HorseshoePanel([0,-1,np.tan(10*np.pi/180)],[0,0,0],1,1,10,10,0,0)
    print('ncontrol = ',H.get_ncontrol())