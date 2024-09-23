#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 16:38:23 2022

@author: leon
"""

import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import os

class IONS:
    """Class for ions containing name, charge, and concentration"""
    def __init__(self, name, charge, concentration):
        # name:String, charge:int, conc:float or array
        self.name = name
        self.charge = charge
        self.conc = concentration
    
    
def calc_debye(ions, er=80):
    e = 1.6 * (10**-19) # C
    KbT = 4.114 * (10**-21) # kg * m^2 * s^-2
    e0 = 8.85418782 * (10**-12) # m^-3 * kg^-1 * s^4 * A^2
    
    sum_charge = 0
    for ion in ions:
        #convert molar concentration to number density
        conc = ion.conc * 6.022*(10**23) * 1000
        sum_charge += (ion.charge * e)**2 * conc
    k = np.sqrt(sum_charge/(e0*er*KbT))
    lambD = k**-1
    
    return k, lambD

def calc_bjerrum(er=80):
    e = 1.6 * (10**-19) # C
    KbT = 4.114 * (10**-21) # kg * m^2 * s^-2
    e0 = 8.85418782 * (10**-12)
    return e**2/(4*np.pi*e0*er*KbT)


def radical_distribution(Universe, ref, sel, dimension="xyz", b=0, e=-1, com=False, relative_rho=True):
    ''' Universe: MDA.universe
        ref: reference point or axis(two points)
        sel: MDA.atomgroups
        b: starting time
        e: ending time
        dimension: "xy" for cylindrical distribution, "xyz" for spherical
        com: calculate distribution using center of mass instead of every atom
        relative_rho: if normalize with respect to bulk concentration
        '''
   
    #calculate the rdf around a point (xyz) or around an axis (xy)
    
    #first determine if com, if true, calculate the COM of each residue
    if com == True:
        ref_com = np.array()
        
    #for 2D, only support vertical Z axis (x=c1, y=c2)
    distri = []
    if dimension == "xyz":
        for ts in Universe.trajectory[b:e]:
            for atom in sel:
                distance = np.linalg.norm(atom.position - ref)
                distri.append(distance)
        if relative_rho == True:
            rho0 = len(sel)/(Universe.dimensions[0]*Universe.dimensions[1]*Universe.dimensions[2]/1000)
        else:
            rho0 = 1/(6.022*10**23)*10**24
        
    if dimension == "xy":
        def lineseg_dist(p, a, b):
            #def a function to calculate the distance from a point to a line
            # normalized tangent vector
            d = np.divide(b - a, np.linalg.norm(b - a))
            # signed parallel distance components
            s = np.dot(a - p, d)
            t = np.dot(p - b, d)
            # clamped parallel distance
            h = np.maximum.reduce([s, t, 0])
            # perpendicular distance component
            c = np.cross(p - a, d)
        
            return np.hypot(h, np.linalg.norm(c)) 
            
        for ts in Universe.trajectory[b:e]:
            for atom in sel:
                distance = lineseg_dist(atom.position, ref[0], ref[1])
                distri.append(distance)
        rho0 = len(sel)/(Universe.dimensions[0]*Universe.dimensions[1]/100)
        if relative_rho == False:
            rho0 = rho0 / (len(sel)/(Universe.dimensions[0]*Universe.dimensions[1]*Universe.dimensions[2]/1000))/(6.022*10**23)*10**24

    hist, edges = np.histogram(
    distri,
    bins=200,
    range=[0,55],
    density=False)
    
    # rho0 = rho* volume of box / area of the cross section of the box
    
    gr = []
    for i in range(len(hist)):
        if hist[i] == 0: 
            gr.append(0)
        else:
            if dimension == 'xy':
                gr.append(float(hist[i])/((2*np.pi*0.0265*i*0.0265)*rho0*len(Universe.trajectory[b:e])))
            else:
                gr.append(float(hist[i])/((4*np.pi*((0.0265*i)**2)*0.0265)*rho0*len(Universe.trajectory[b:e])))
                
    plt.plot(edges[:-1], gr)
    
    return gr, edges



