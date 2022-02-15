"""
Please contact Gianluca Meneghello <gianmail@gmail.com> for further information

Copyright (C) 2020 Gianluca Meneghello 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import sys
from numpy import *
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.interpolate import InterpolatedUnivariateSpline 
from scipy.integrate import ode

def readDict(filename):
  ''' reads a file and return a dictionary of name-value combinations
  anything line not starting with a letter is discarded'''
  d = {}
  with open(filename) as f:
    for line in f:
      if (line[0].isalpha()):
        splitvec = line.strip().split(None,2)
        d[splitvec[0]] = float(splitvec[1])
  return d

##########################
### HYDRODYNAMIC SETUP ###
##########################
def RrhSetUp(geom,phys):
  '''Hull residuary resistance'''
  FnD = array(r_[.10:.61:.05])

  coeff = array( [[-0.0014 ,  0.0403 ,  0.0470 , -0.0227 , -0.0119 ,  0.0061 , -0.0086 , -0.0307 , -0.0553] , 
                  [ 0.0004 , -0.1808 ,  0.1793 , -0.0004 ,  0.0097 ,  0.0118 , -0.0055 ,  0.1721 , -0.1728] , 
                  [ 0.0014 , -0.1071 ,  0.0637 ,  0.0090 ,  0.0153 ,  0.0011 ,  0.0012 ,  0.1021 , -0.0648] , 
                  [ 0.0027 ,  0.0463 , -0.1263 ,  0.0150 ,  0.0274 , -0.0299 ,  0.0110 , -0.0595 ,  0.1220] , 
                  [ 0.0056 , -0.8005 ,  0.4891 ,  0.0269 ,  0.0519 , -0.0313 ,  0.0292 ,  0.7314 , -0.3619] , 
                  [ 0.0032 , -0.1011 , -0.0813 , -0.0382 ,  0.0320 , -0.1481 ,  0.0837 ,  0.0223 ,  0.1587] , 
                  [-0.0064 ,  2.3095 , -1.5152 ,  0.0751 , -0.0858 , -0.5349 ,  0.1715 , -2.4550 ,  1.1865] , 
                  [-0.0171 ,  3.4017 , -1.9862 ,  0.3242 , -0.1450 , -0.8043 ,  0.2952 , -3.5284 ,  1.3575] , 
                  [-0.0201 ,  7.1576 , -6.3304 ,  0.5829 ,  0.1630 , -0.3966 ,  0.5023 , -7.1579 ,  5.2534] , 
                  [ 0.0495 ,  1.5618 , -6.0661 ,  0.8641 ,  1.1702 ,  1.7610 ,  0.9176 , -2.1191 ,  5.4281] , 
                  [ 0.0808 , -5.3233 , -1.1513 ,  0.9663 ,  1.6084 ,  2.7459 ,  0.8491 ,  4.7129 ,  1.1089]])

  vect = array( [geom['LWL']/geom['DIVCAN']**(1./3.) , geom['XFB']/geom['LWL'] , geom['CPL'] , geom['DIVCAN']**(2./3.)/geom['AW'] , geom['BWL']/geom['LWL'] , geom['DIVCAN']**(2./3.)/geom['SC'] , geom['XFB']/geom['XFF'] , (geom['XFB']/geom['LWL'])**2. , geom['CPL']**2.])

  RrhD = (geom['DIVCAN']*phys['g']*phys['rho_w']) * dot( coeff , vect ) * (geom['DIVCAN']**(1./3.)/geom['LWL'])

  Rrh = InterpolatedUnivariateSpline(FnD,RrhD) 

  return Rrh

def RrhHSetUp(geom,phys):
  '''Change in hull residuary resistance due to heel'''
  FnD = array(r_[.10:.61:.05])

  coeff = array( [[  0      ,  0      ,  0      ,  0      ,  0      ,  0     ] ,
                  [  0      ,  0      ,  0      ,  0      ,  0      ,  0     ] ,
                  [  0      ,  0      ,  0      ,  0      ,  0      ,  0     ] ,
                  [ -0.0268 , -0.0014 , -0.0057 ,  0.0016 , -0.0070 , -0.0017] ,
                  [  0.6628 , -0.0632 , -0.0699 ,  0.0069 ,  0.0459 , -0.0004] ,
                  [  1.6433 , -0.2144 , -0.1640 ,  0.0199 , -0.0540 , -0.0268] ,
                  [ -0.8659 , -0.0354 ,  0.2226 ,  0.0188 , -0.5800 , -0.1133] ,
                  [ -3.2715 ,  0.1372 ,  0.5547 ,  0.0268 , -1.0064 , -0.2026] ,
                  [ -0.1976 , -0.1480 , -0.6593 ,  0.1862 , -0.7489 , -0.1648] ,
                  [  1.5873 , -0.3749 , -0.7105 ,  0.2146 , -0.4818 , -0.1174] ,
                  [  1.5873 , -0.3749 , -0.7105 ,  0.2146 , -0.4818 , -0.1174]] ) / 1000.


  vect = array( [ 1 , geom['LWL']/geom['BWL'] , geom['BWL']/geom['TCAN'] , (geom['BWL']/geom['TCAN'])**2 , geom['XFB']/geom['LWL'] , (geom['XFB']/geom['LWL'])**2] ) 

  RrhH20D =  (geom['DIVCAN']*phys['g']*phys['rho_w']) * dot( coeff , vect )

  RrhH20 = InterpolatedUnivariateSpline(FnD,RrhH20D) 

  def RrhH(Fn,phi): return RrhH20(Fn) * 6. * abs(phi)**1.7

  return RrhH

def RrkSetUp(geom,phys):
  '''Keel viscous resistance'''
  FnD = array(r_[.10:.61:.05])
  
  coeff = array( [[  0       , 0       ,  0       ,  0      ],
                  [  0       , 0       ,  0       ,  0      ],
                  [ -0.00104 , 0.00172 ,  0.00117 , -0.00008],
                  [ -0.00550 , 0.00597 ,  0.00390 , -0.00009],
                  [ -0.01110 , 0.01421 ,  0.00069 ,  0.00021],
                  [ -0.00713 , 0.02632 , -0.00232 ,  0.00039],
                  [ -0.03581 , 0.08649 ,  0.00999 ,  0.00017],
                  [ -0.00470 , 0.11592 , -0.00064 ,  0.00035],
                  [  0.00553 , 0.07371 ,  0.05991 , -0.00114],
                  [  0.04822 , 0.00660 ,  0.07048 , -0.00035],
                  [  0.01021 , 0.14173 ,  0.06409 , -0.00192]] )

  vect = array ( [ 1 , geom['T']/geom['BWL'] , (geom['T']-geom['ZCBK'])/geom['DVK']**(1./3.) , geom['DIVCAN']/geom['DVK'] ] ) 

  RrkD = (geom['DVK']*phys['g']*phys['rho_w']) * dot( coeff , vect )

  Rrk = InterpolatedUnivariateSpline(FnD,RrkD) 

  return Rrk

def RrkHSetUp(geom,phys):
  '''Change in keel resiudary resistance due to heel'''
  coeff = array( [ -3.5837 , -0.0518 , 0.5958 , 0.2055] )
  vect  = array( [geom['TCAN']/geom['T'] , geom['BWL']/geom['TCAN'] , (geom['BWL']*geom['TCAN'])/(geom['T']*geom['TCAN']) , geom['LWL']/geom['DIVCAN']**1/3 ] )
  Ch = dot( coeff , vect )
 
  def RrkH(Fn,phi):  return geom['DVK']*phys['rho_w']*phys['g']*Ch*Fn**2*abs(phi)
 
  return RrkH

def RiSetUp(geom,phys):
  '''Keel induced resistance'''
  coeffA = array([[3.7455 , -3.6246 , 0.0589 , -0.0296],
                  [4.4892 , -4.8454 , 0.0294 , -0.0176],
                  [3.9592 , -3.9804 , 0.0283 , -0.0075],
                  [3.4891 , -2.9577 , 0.0250 , -0.0272]] )

  vectA = array([ geom['TCAN']/geom['T'] , (geom['TCAN']/geom['T'])**2 , geom['BWL']/geom['TCAN'] , geom['CHTPK']/geom['CHRTK']] )
  
  Tegeo = dot(coeffA,vectA)

  def Ri(Fn,phi,Fside): 
    phiD = array( [0., 10., 20., 30.] ) * pi / 180. 

    coeffB = array([[1.2306 , -0.7256],
                    [1.4231 , -1.2971],
                    [1.545  , -1.5622],
                    [1.4744 , -1.3499]] )

    vectB = array([1., Fn]);

    TeFn  = dot(coeffB,vectB)

    TeD = geom['T'] *  Tegeo * TeFn;

    Te = InterpolatedUnivariateSpline(phiD,TeD);

    # calculating Ri, the induced resistance...
    Fheel = Fside/cos(phi)

    V = Fn*sqrt(phys['g']*geom['LWL'])

    return Fheel**2 / (pi * Te(abs(phi))**2 * 0.5 * phys['rho_w'] * V**2)

  return Ri

def RvhSetUp(geom,phys):
  '''Hull viscous resistance'''
  def Rvh(Fn):
    V = Fn*sqrt(phys['g']*geom['LWL'])
    Rn = geom['LWL'] * 0.7 * V / phys['ni_w']
    Cf = 0.075 / ( log10(Rn) - 2 )**2
    Rfh = 1./2. * phys['rho_w'] * V**2 * geom['SC'] * Cf
    return Rfh * geom['HULLFF']

  return Rvh

def RvhHSetUp(geom,phys):
  '''Change in hull viscous resistance due to heel'''
  phiD = array( [0. , 5. , 10. , 15. , 20. , 25. , 30. , 35.] ) * pi/180.

  coeff = array( [[  0.000 ,  0.000 ,  0.000 , 0.000],
                  [ -4.112 ,  0.054 , -0.027 , 6.329],
                  [ -4.522 , -0.132 , -0.077 , 8.738],
                  [ -3.291 , -0.389 , -0.118 , 8.949],
                  [  1.850 , -1.200 , -0.109 , 5.364],
                  [  6.510 , -2.305 , -0.066 , 3.443],
                  [ 12.334 , -3.911 ,  0.024 , 1.767],
                  [ 14.648 , -5.182 ,  0.102 , 3.497]] )

  vect = array( [1. , geom['BWL']/geom['TCAN'] , (geom['BWL']/geom['TCAN'])**2. , geom['CMS'] ] ) 

  SCphiD = geom['SC'] * ( 1. + 1./100.*dot(coeff,vect))
 
  SCphi = InterpolatedUnivariateSpline(phiD,SCphiD)

  def RvhH(Fn,phi):
    V = Fn*sqrt(phys['g']*geom['LWL'])
    Rn = geom['LWL'] * 0.7 * V / phys['ni_w'];
    Cf = 0.075 / ( log10(Rn) - 2. )**2.
    RfhH = 1./2. * phys['rho_w'] * V**2 * Cf * (SCphi(abs(phi)) - geom['SC'])
    return RfhH * geom['HULLFF']  #does it make sense to use the same hull form factor both for the upright and the heeled hull?

  return RvhH

def RvkSetUp(geom,phys):
  '''Keel viscous resistance'''
  def Rvk(Fn):
    V = Fn*sqrt(phys['g']*geom['LWL'])
    Rn = geom['CHMEK'] * V / phys['ni_w'];
    Cf = 0.075 / ( log10(Rn) - 2 )**2;
    Rfk = 1./2. * phys['rho_w'] * V**2 * geom['SK'] * Cf;
    return Rfk * geom['KEELFF']

  return Rvk

def RvrSetUp(geom,phys):
  '''Rudder viscous resistance'''
  def Rvr(Fn):
    V = Fn*sqrt(phys['g']*geom['LWL'])
    Rn = geom['CHMER'] * V / phys['ni_w'];
    Cf = 0.075 / ( log10(Rn) - 2 )**2;
    Rfr = 1./2. * phys['rho_w'] * V**2 * geom['SR'] * Cf;
    return Rfr * geom['RUDDFF'];

  return Rvr

#########################
### AERODYNAMIC SETUP ###
#########################
def aeroSetUp(geom,phys):
  '''calculation of sails' area and ZCE (above waterline)'''
  AM = 0.5 * geom['P'] * geom['E'] *geom['MROACH']             #mainsail area
  AJ = 0.5 * sqrt( geom['I']**2 + geom['J']**2) * geom['LPG']   #jib area
  AS = 1.15 * geom['SL'] * geom['J']                         #spi area
  AF = 0.5 * geom['I'] * geom['J']                           #foretriangle area
  
  ZCEM = 0.39 * geom['P'] + geom['BAD']
  ZCEJ = 0.39 * geom['I']
  ZCES = 0.59 * geom['I']

  # calculating Cl and Cd
  
  # LE_Cl-LE_Cd -> Hazer Cl-Cd coefficients, 1999
  # [Apparent wind angle, Main, Jib, spi]
  if (geom['MFLB'] == 1):
    # Revision for full-length battens in main; added zeros at alfa = 0
    LE_Cl = array( [[0.   , 0.    , 0.  , 0.  ],
                    [27.  , 1.725 , 1.5 , 0.  ],
                    [50.  , 1.5   , 0.5 , 1.5 ],
                    [80.  , 0.95  , 0.3 , 1.0 ],
                    [100. , 0.85  , 0.0 , 0.85],
                    [180. , 0.    , 0.  , 0.  ]] )
  else:
    LE_Cl = array( [[ 0.  , 0.   , 0.  , 0.   ],
                    [ 27. , 1.5  , 1.5 , 0.   ],
                    [ 50. , 1.5  , 0.5 , 1.5  ],
                    [ 80. , 0.95 , 0.3 , 1.0  ],
                    [ 100., 0.85 , 0.0 , 0.85 ],
                    [ 180., 0.   , 0.  , 0.   ]] )
  
  LE_Cdp = array( [[ 0.  , 0.   , 0.   , 0.   ],
                   [ 27. , 0.02 , 0.02 , 0.   ],
                   [ 50. , 0.15 , 0.25 , 0.25 ],
                   [ 80. , 0.8  , 0.15 , 0.9  ],
                   [ 100., 1.0  , 0.0  , 1.2  ],
                   [ 180., 0.9  , 0.0  , 0.66 ]] )

  # Changed from spline to pchip to avoid oscillations that produced negative
  # drag and other non-physical behavior.
  Cl_M = InterpolatedUnivariateSpline(LE_Cl[:,0]*pi/180.,LE_Cl[:,1])
  Cl_J = InterpolatedUnivariateSpline(LE_Cl[:,0]*pi/180.,LE_Cl[:,2])
  Cl_S = InterpolatedUnivariateSpline(LE_Cl[:,0]*pi/180.,LE_Cl[:,3])
  
  Cdp_M = InterpolatedUnivariateSpline(LE_Cdp[:,0]*pi/180.,LE_Cdp[:,1])
  Cdp_J = InterpolatedUnivariateSpline(LE_Cdp[:,0]*pi/180.,LE_Cdp[:,2])
  Cdp_S = InterpolatedUnivariateSpline(LE_Cdp[:,0]*pi/180.,LE_Cdp[:,3])

  def aero(V,phi,F,V_tw,alfa_tw):
    
    # rotating the wind vector in the mast's perpendicular plane... this is vpp_aero in the matlab version
    V1 = V + V_tw * cos(alfa_tw);
    V2 = V_tw * abs(sin(alfa_tw)) * cos(phi);
    
    V_eff = sqrt( V1**2 + V2**2 );
    alfa_eff = arctan2(V2,V1);
    
    if (geom['SAILSET'] == 1):   # main only
      AN = AM
      Cl = Cl_M(alfa_eff)
      Cdp = Cdp_M(alfa_eff)
      ZCE = ZCEM
    elif (geom['SAILSET'] == 3):  # main and jib
      AN = AF + AM
      Cl = ( Cl_M(alfa_eff) * AM + Cl_J(alfa_eff)* AJ ) / AN
      Cdp = ( Cdp_M(alfa_eff) * AM + Cdp_J(alfa_eff) * AJ ) / AN
      ZCE = (ZCEM * alfa_eff * AM + ZCEJ * AJ) / (AM + AJ)
    elif  (geom['SAILSET'] == 5): # main and spi
      AN = AF + AM
      Cl = ( Cl_M(alfa_eff) * AM + Cl_S(alfa_eff) * AS ) / AN;
      Cdp = ( Cdp_M(alfa_eff) * AM + Cdp_S(alfa_eff) * AS ) / AN;
      ZCE = (ZCEM * AM + ZCES * AS) / (AM + AS);
    elif (geom['SAILSET'] == 7):  # main and jib and spi
      AN = AF + AM
      Cl = ( Cl_M(alfa_eff) * AM + Cl_J(alfa_eff) * AJ + Cl_S(alfa_eff) * AS ) / AN;
      Cdp = ( Cdp_M(alfa_eff) * AM + Cl_J(alfa_eff) * AJ + Cdp_S(alfa_eff) * AS ) / AN;
      ZCE = (ZCEM * AM + ZCEJ * AJ + ZCES * AS) / (AM + AJ + AS);
    else:
      print('VALID VALUES FOR SAILSET ARE: 1, 3, 5, 7')

    Cl = Cl*F

    if (alfa_eff < 45*pi/180.):
      AR = (1.1*( geom['EHM'] + geom['AVGFREB']))**2/AN;
    else:
      AR = (1.1*( geom['EHM'] ))**2/AN;

    CdI = Cl**2 * ( 1/(pi*AR) + 0.005); #induced resistance
    Cd0 = 1.13  * ( (geom['B']*geom['AVGFREB']) + (geom['EHM']*geom['EMDC']) ) / AN;
    
    Cd = Cdp + Cd0 + CdI;

    L = 0.5 * phys['rho_a'] * V_eff**2 * AN * Cl
    D = 0.5 * phys['rho_a'] * V_eff**2 * AN * Cd
    
    Fdrive = L * sin(alfa_eff) - D * cos(alfa_eff)
    Fheel = L * cos(alfa_eff) + D * sin(alfa_eff)
    #if Fdrive < 0; Fdrive = 0; end
    #if Fheel < 0; Fheel = 0; end
    Mheel = Fheel*(ZCE + geom['T'] - geom['ZCBK']) # attenzione: il centro di spinta della deriva e' stato messo nel centro di galleggiamento. l'ipotesi e' corretta?
    Fside = Fheel*cos(phi)

    return (Fdrive,Fside,Mheel)
  
  return aero

#####################
### THE VPP CLASS ###
#####################
class vpp:
  "The vpp class"
  
  V2Fn = lambda self,V  : V/sqrt(self.phys['g']*self.geom['LWL'])
  Fn2V = lambda self,Fn : Fn*sqrt(self.phys['g']*self.geom['LWL'])

  def __init__(self,boatname):

    self.geom = readDict(boatname+'.geom')     #loading dictionaryfrom files
    self.sailset = {
        1: 'main only',
        3: 'main and jib',
        5: 'main and spinnaker',
        7: 'main jib and spinnaker'
        }[self.geom['SAILSET']]
    self.phys = {
        'rho_w': 1025.9,        # water density
        'ni_w': 1.18838E-6,     # water viscosity
        'rho_a': 1.125,         # air density
        'g': 9.80665            # gravitational accelerations
        }
    
    # hydrodynamic initialization, all angles are in radians 
    self.Rrh  = RrhSetUp  ( self.geom , self.phys ) 
    self.RrhH = RrhHSetUp ( self.geom , self.phys ) 
    self.Rrk  = RrkSetUp  ( self.geom , self.phys ) 
    self.RrkH = RrkHSetUp ( self.geom , self.phys ) 
    self.Ri   = RiSetUp   ( self.geom , self.phys ) 
    self.Rvh  = RvhSetUp  ( self.geom , self.phys ) 
    self.RvhH = RvhHSetUp ( self.geom , self.phys ) 
    self.Rvk  = RvkSetUp  ( self.geom , self.phys ) 
    self.Rvr  = RvrSetUp  ( self.geom , self.phys ) 

    # aerodynamic inistialization
    self.aero = aeroSetUp ( self.geom , self.phys ) 

  def Rtot(self,Fn,phi,Fside): 
    return self.Ri(Fn,phi,Fside) + \
           self.Rrh(Fn)          + \
           self.RrhH(Fn,phi)     + \
           self.Rrk(Fn)          + \
           self.RrkH(Fn,phi)     + \
           self.Rvh(Fn)          + \
           self.RvhH(Fn,phi)     + \
           self.Rvk(Fn)          + \
           self.Rvr(Fn)

  def Mright(self,phi,b):
    M1 = (self.geom['KM'] - self.geom['KG']) * sin(phi) * self.phys['g'] * self.phys['rho_w'] * (self.geom['DIVCAN'] + self.geom['DVK'])
    M2 = self.geom['MMVBLCRW'] * self.phys['g'] * b * cos(phi)
    Mright = M1+M2

    return Mright

  def force(self,V,phi,b,F,V_tw,alfa_tw):
    (Fdrive,Fside,Mheel) = self.aero(V,phi,F,V_tw,alfa_tw)
    Fn = self.V2Fn(V) #V/sqrt(self.phys['g']*self.geom['LWL'])
    return Fdrive - self.Rtot(Fn,phi,Fside)

  def momentum(self,V,phi,b,F,V_tw,alfa_tw):
    (Fdrive,Fside,Mheel) = self.aero(V,phi,F,V_tw,alfa_tw)
    return Mheel - self.Mright(phi,b)

  #############################
  ## THE ACTUAL OPTIMIZATION ##
  #############################
  def solve(self,V_tw,alfa_tw,x0 = array([2.,0.,1.,1.]) ):
    cons = (
        {'type': 'eq', 'fun' : lambda x: self.force   (x[0],x[1],x[2],x[3],V_tw,alfa_tw) },   # constrain on force
        {'type': 'eq', 'fun' : lambda x: self.momentum(x[0],x[1],x[2],x[3],V_tw,alfa_tw) }    # constrain on momentum
        )

    bnds = [
        (self.Fn2V(0.1), self.Fn2V(0.6)),   # bounds on velocity (m/s)
        ( 0, 50*pi/180),                    # bounds on heeling angle (radians)
        (0, 1),                            # bounds on crew positions (m)
        (.4, 1)                             # bounds on sail flattening parameter
        ]

    # maximize the squared of the velocity (i.e. minimize the negative velocity)
    res = minimize(lambda x: -x[0]**2, x0 ,constraints=cons, bounds=bnds, method='SLSQP', options={'disp': False})

    return res

  #########################
  ## THE DYNAMICAL MODEL ##
  #########################
  def race(self):
    y0, t0 = array([0.,3,40.*pi/180,0.]), 0

    def f(t, y, arg1):

      F = 1.
      V_tw = 10.
      alfa_tw = 90.*pi/180
      b=3
      mass = 1150 + 6*80
      massI = mass*8
      r = 1 # damping parameters for angular momentum equation

      x, V, phi, phidot = y
      (Fdrive,Fside,Mheel) = self.aero(V,phi,F,V_tw,alfa_tw)
      Fn = self.V2Fn(V)

      res = array([ V , ( Fdrive - self.Rtot(Fn,phi,Fside) )/mass, phidot-r*phi ,  (Mheel - self.Mright(phi,b) )/massI ])

      return res

    r = ode(f).set_integrator('lsoda')
    r.set_initial_value(y0, t0).set_f_params(2.0).set_jac_params(2.0)
    t1 = 500
    dt = .1
    plt.close('run')
    fig,ax1 = plt.subplots(1,1,num='run')
    ax2 = ax1.twinx()
    ax2.grid()
    while r.successful() and r.t < t1:
      sol = r.integrate(r.t+dt)
      print(r.t+dt, sol)
      ax1.plot(r.t+dt, sol[1],'.k')
      ax2.plot(r.t+dt, sol[2]*180/pi,'.r')
      plt.pause(.01)

#####################
### MAIN FUNCTION ###
#####################
if __name__ == "__main__":

  # read command line arguments
  try :
    boatname = sys.argv[1]
    V_tws = [ float(tmp) for tmp in sys.argv[2:] ]
  except IndexError as idx:
    raise RuntimeError('Missing input argument. To call the vpp use:\n\npython vpp.py <boatname> <windSpeed> [<windSpeed> ...]\n\n') from idx

  # initialize class
  boat = vpp(boatname)

  # initialize figure and output
  fig,axs = plt.subplots(1,4,figsize=(12,4),subplot_kw=dict(projection="polar"))
  fig.suptitle('boat name: '+boatname+'\nsailset: '+boat.sailset)

  print("# Computing polar for ",boatname,' and wind speed ',V_tws,'[m/s]')
  print("%12s %12s %12s %12s %12s %12s | %18s %18s %30s"%('V_tw','windAngle','boatSpeed','heel','crewArm','sailParam','residualForce','residualMomentum','solverstatus'))
  print("%12s %12s %12s %12s %12s %12s | %18s %18s %30s"%('[m/s]','[degrees]','[m/s]','[degress]','[m]','[-]','[N]','[N]','[-]'))
  
  # solve for all wind speed
  for V_tw in V_tws:
    res = {}
    x0 = array([ 4.,0.,0.,4.])    # V, phi, b, F
    polar = []
    windAngles = r_[ 25:181:5]*pi/180
    # solve for all wind angles
    for i,alfa_tw in  enumerate(windAngles):
      res = boat.solve(V_tw,alfa_tw,x0)
      polar.append(res.x.copy())
      x0 = res.x

      # verify that equilibrium has been reached
      (Fdrive,Fside,Mheel) =  boat.aero(res.x[0],res.x[1],res.x[3],V_tw,alfa_tw)
      Rtot = boat.Rtot(boat.V2Fn(res.x[0]),res.x[1],Fside)
      Mright = boat.Mright(res.x[1],res.x[2])
      
      print("%12.0f %12.0f %12.1f %12.1f %12.1f %12.1f | %+18.2e %+18.2e %30s"%(V_tw,degrees(alfa_tw),res.x[0],degrees(res.x[1]),res.x[2],res.x[3],(Fdrive-Rtot)/Rtot,(Mright-Mheel)/Mheel,res.message))

    #### PLOTTING ###
    polar = asarray(polar)
    ax= axs[0]
    ax.set_title('boat speed [m/s]')
    ax.set_theta_zero_location("N")
    ax.plot(windAngles,polar[:,0],label=str(V_tw))
    ax.set_xticklabels([])

    ax = axs[1]
    ax.set_title('boat heel [degress]')
    ax.set_theta_zero_location("N")
    ax.plot(windAngles,degrees(polar[:,1]))
    ax.set_xticklabels([])

    ax = axs[2]
    ax.set_title('crew arm [m]')
    ax.set_theta_zero_location("N")
    ax.plot(windAngles,polar[:,2])
    ax.set_xticklabels([])

    ax = axs[3]
    ax.set_title('sail trim [m]')
    ax.set_theta_zero_location("N")
    ax.plot(windAngles,polar[:,3])
    ax.set_xticklabels([])

  fig.legend(title='Wind speed [m/s]',loc=8,ncol=len(V_tws))
  plt.tight_layout()
  figname = boatname+'.png'
  plt.savefig(figname,bbox_inches='tight')
  print('figures saved in %s'%figname)
  print('close all windows to exit')
  plt.show()
