#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class for modeling Unified Engineering Flight Competition wings using a vortex 
lattice method

Author: David Darmofal, MIT
Date: April 2, 2021
"""

import numpy as np
import HorseshoePanel as HP
import matplotlib.pyplot as plt

class UEFC_wing(object):
    def __init__(self, b=None, croot=None, ctip=None, agroot=None, agtip=None, dihedral=None):
 
        self.b        = b
        self.croot    = croot
        self.ctip     = ctip
        self.agroot   = agroot
        self.agtip    = agtip
        self.dihedral = dihedral
    
        Nsemi = 50
        
        al0 = -6.0 # This is a reasonable value for all UEFC airfoils
        pi = np.pi
        cmax = max(croot,ctip)
        self.HPlist = []
        for i in range(2*Nsemi):
            tha = 0.5*i*pi/Nsemi
            thb = tha + 0.5*pi/Nsemi
            fa = -np.cos(tha)
            fb = -np.cos(thb)
            ya  = 0.5*b*fa
            yb  = 0.5*b*fb
            afa = np.abs(fa)
            afb = np.abs(fb)
            ca  = croot + (ctip-croot)*afa
            cb  = croot + (ctip-croot)*afb
            aga = agroot + (agtip-agroot)*afa
            agb = agroot + (agtip-agroot)*afb
            xa  = cmax-ca
            xb  = cmax-cb
            za  = 0.5*b*afa*np.tan(dihedral*pi/180)
            zb  = 0.5*b*afb*np.tan(dihedral*pi/180)
            Xa  = np.array([xa, ya, za])
            Xb  = np.array([xb, yb, zb])
            pan = HP.HorseshoePanel(Xa, Xb, ca, cb, aga, agb, al0, al0)
            self.HPlist.append(pan)
        
    
    def plotgeom(self):
        fig, axs = plt.subplots(3,1,sharex=True)
        axs_p = axs[0] # planform plot
        axs_f = axs[1] # front plot
        axs_g = axs[2] # geometric twist
    
        for pan in self.HPlist:
            XLEa = pan.get_XLEa()
            XLEb = pan.get_XLEb()
            ca   = pan.get_ca()
            cb   = pan.get_cb()
            aga  = pan.get_aga()
            agb  = pan.get_agb()
            
            axs_p.plot([XLEa[1],XLEb[1],XLEb[1],XLEa[1],XLEa[1]],[XLEa[0],XLEb[0],XLEb[0]+cb,XLEa[0]+ca,XLEa[0]],color='black')
            
            axs_f.plot([XLEa[1],XLEb[1]],[XLEa[2],XLEb[2]],color='black')
            
            axs_g.plot([XLEa[1],XLEb[1]],[aga, agb],color='black')
        
        axs_p.invert_yaxis()    
        axs_p.axis('equal')
        axs_p.set_ylabel('$x$')
        axs_f.axis('equal')
        axs_f.set_ylabel('$z$')
        axs_g.set_ylabel(r'$\alpha_g$')
        axs_g.set_xlabel('$y$')
        axs_p.set_title('$S =$ {:.3f}, $AR =$ {:.2f}'.format(self.get_S(),self.get_AR()))
        
        return (fig, axs)
        
    
    def get_S(self):
        return 0.5*(self.croot+self.ctip)*self.b
    
    def get_AR(self):
        return (self.b**2)/self.get_S()
    
    def solve(self, alpha=None, CL=None):
        
        pi = np.pi
        if not (alpha is None):
            alphar = alpha*pi/180.
            alphad = alpha
            Vinf = np.array([1.0, 0., alphar]) # small angle approimation
        elif (CL is None):
            print('Error: you must define alpha or CL')
        
        N = len(self.HPlist)
        A = np.zeros((N,N))
        b = np.zeros(N)
        for i in range(N):
            pani = self.HPlist[i]
            nci  = pani.get_ncontrol()
            Xci  = pani.get_Xcontrol()
            
            for j in range(N):
                panj = self.HPlist[j]
                A[i,j] = panj.get_unitV(Xci) @ nci
        
        if not (alpha is None):
            for i in range(N):
                pani = self.HPlist[i]
                nci  = pani.get_ncontrol()
                Xci  = pani.get_Xcontrol()
                b[i] = -np.matmul(Vinf,nci)
            
            G = np.linalg.solve(A,b)
        
        else:
            Vinf = np.array([1.0, 0.0, 0.0])
            for i in range(N):
                pani = self.HPlist[i]
                nci  = pani.get_ncontrol()
                Xci  = pani.get_Xcontrol()
                b[i] = -np.matmul(Vinf,nci)           
                
            G0 = np.linalg.solve(A,b)
            CL0, CD0, e0, clmax0 = self.calc_aeroperf(G0)
            
            Vinf = np.array([1.0, 0.0, 0.1])
            for i in range(N):
                pani = self.HPlist[i]
                nci  = pani.get_ncontrol()
                Xci  = pani.get_Xcontrol()
                b[i] = -np.matmul(Vinf,nci) 
                
            G1 = np.linalg.solve(A,b)
            CL1, CD1, e1, clmax1 = self.calc_aeroperf(G1)

            G = G0 + (CL-CL0)/(CL1-CL0)*(G1-G0)
            alphad = 0.1*180./pi*(CL-CL0)/(CL1-CL0)
        
        return G, alphad
    
    
    def plotcl(self, G, plotclccbar=False):
    
        cl, y = self.calccldist(G)
        fig, axs = plt.subplots(1,1,sharex=True)
        axs.plot(2*y/self.b,cl,'r-',label='$c_l$')
        cbar = 0.5*(self.croot + self.ctip)
        clccbar = 2*G/cbar
        if plotclccbar:
            axs.plot(2*y/self.b,clccbar,'b-',label=f'$c_l\\, c/\\overline{{c}}$')
        
        clmax = max(cl)
        CL, CDi, e0, clmax0 = self.calc_aeroperf(G)
        axs.set_title('$C_L =$ {:.2f}, $\\mathrm{{max}}(c_l) =$ {:.2f}, $C_{{D_i}} =$ {:.4f}, $e_0 =$ {:.2f}'.format(CL,clmax,CDi,e0))
        axs.set_xlabel('$y/(b/2)$')
        axs.grid(True)
        axs.legend()
        
        return (fig, axs)
    
    
    def calccldist(self, G):
        N = len(self.HPlist)
        y = np.zeros(N)
        c = np.zeros(N)
        for i in range(N):
            pan = self.HPlist[i]
            XLEa = pan.get_XLEa()
            XLEb = pan.get_XLEb()
            ca   = pan.get_ca()
            cb   = pan.get_cb()
            
            y[i] = 0.5*(XLEa[1] + XLEb[1])
            c[i] = 0.5*(ca+cb)
        
        cl = 2*G/c
        
        return cl, y
    
    
    def calc_aeroperf(self, G):
        
        # Calculate maximum sectional cl
        cl, y = self.calccldist(G)
        clmax = max(cl)
        
        # Calculate lift
        N = len(self.HPlist)
        L = 0
        for i in range(N):
            pan = self.HPlist[i]
            Xa = pan.get_Xbounda()
            Xb = pan.get_Xboundb()
            dy = Xb[1]-Xa[1]
            L = L + G[i]*dy

        # Calculate induced drag
        Di = 0.
        for j in range(N+1):
            if (j==0):
                dG = G[0]
            elif (j==N):
                dG = -G[-1]
            else:
                dG = G[j]-G[j-1]
            
            if (j==N):
                Xj = self.HPlist[N-1].get_Xboundb()
            else:
                Xj = self.HPlist[j].get_Xbounda()
                
            yj = Xj[1]
            zj = Xj[2]
            
            pi = np.pi
            for i in range(N):
                pan = self.HPlist[i]
                Xa = pan.get_Xbounda()
                Xb = pan.get_Xboundb()
                dy = Xb[1]-Xa[1]
                dz = Xb[2]-Xa[2]
                ds = np.sqrt(dy**2 + dz**2)
                ny = -dz/ds
                nz =  dy/ds
                Xci = pan.get_Xcontrol()
                yi = Xci[1]
                zi = Xci[2]
                dyij = yi - yj
                dzij = zi - zj
                dij2 = dyij**2 + dzij**2
                Vn = dG/(2*pi)*(dzij*ny - dyij*nz)/dij2                
                Di = Di + Vn*G[i]*ds
                
                
        Di = -0.5*Di
            
        S = self.get_S()    
        CL = 2*L/S
        CDi = 2*Di/S
        
        e0 = CL**2/(pi*self.get_AR()*CDi)
        return CL, CDi, e0, clmax
            
    
if __name__ == "__main__":
    PV = UEFC_wing(1.5, 0.2, 0.1, 0., -5.0, 10.)
    PV.plotgeom()
    G, alphad = PV.solve(CL=0.85)
    PV.plotcl(G)
    CL, CDi, e0, clmax = PV.calc_aeroperf(G)
    print('CL     = {:.2f}'.format(CL))
    print('clmax  = {:.2f}'.format(clmax))
    print('CDi    = {:.4f}'.format(CDi))
    print('e0     = {:.2f}'.format(e0))



    
    
    
