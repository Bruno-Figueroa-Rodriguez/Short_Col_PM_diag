# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


from units import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#Defining Local geometry

b = 16*inch
h = 12*inch
CC_to_CL_long = 2.5*inch

#Defining Material Properties

fpc = 4000*psi
Ec = 57000*(fpc/psi)**0.5*psi
εcu = 0.003
ß1 = 0.85-(0.05*(fpc-4000*psi))/1000

fy = 60*ksi
Es = 29000*ksi



#Defining Reinforcement


num_lay = 2
bar_per_lay = 3
bar_size = '#8'



def create_rebar_array(num_lay,bar_per_lay,bar_size,CC_to_CL_long):
    rebar_array = []
    for i in range(num_lay):
        rebar_array.append([bar_per_lay,bar_size,CC_to_CL_long+i*(h-2*CC_to_CL_long)/(num_lay-1)])
    
    return rebar_array


def find_strain(Z,rebar_array,εcu):
  εs1 = Z*fy/Es
  c = c = εcu*((h-CC_to_CL_long)/(εcu-εs1))
  a = c*ß1
  if a > h:
    a = h
    
  #print(εs1)
  
  if εs1 <= -εcu-fy/Es:
      φ = 0.9
  elif εs1 < -fy/Es and εs1 > -fy/Es-εcu:
      φ = 0.65+0.25*(-εs1-fy/Es)/εcu
  elif εs1 < 0 and εs1 >= -fy/Es:
      φ = 0.65
  else:
      φ = 0.65


  strains = []

  for reb in rebar_array:

    strains.append(εcu*(c-reb[2])/c)

  return strains,a,φ



###
def find_rebar_force(strains,rebar_array,fpc=4*ksi,fy=60*ksi,Es=29000*ksi):
  forces = []
  #print(len(strains))
  for strain_count in range(len(strains)):

    if strains[strain_count] >= fy/Es:
      forces.append((fy-0.85*fpc)*rebar_array[strain_count][0]*rebars[rebar_array[strain_count][1]]['area'])

    elif strains[strain_count] <= -fy/Es:
      forces.append((-fy)*rebar_array[strain_count][0]*rebars[rebar_array[strain_count][1]]['area'])

    elif strains[strain_count] >0:
      forces.append(((strains[strain_count]*Es)-0.85*fpc)*rebar_array[strain_count][0]*rebars[rebar_array[strain_count][1]]['area'])

    else:
      forces.append((strains[strain_count]*Es*rebar_array[strain_count][0])*rebars[rebar_array[strain_count][1]]['area'])

  return forces
    


def find_PM(a,b,h,forces, rebar_array,φ, fpc=4*ksi):
    """
    This generates a singular point for PM diagram in KIP, KIP-FT
    """
  
    C_c = 0.85*fpc*a*b #Compressive force Kips

    Pn = sum(forces)+C_c
    
    Mn = []
    for force_count in range(len(forces)):
        Mn.append(forces[force_count]*(h/2-rebar_array[force_count][2]))
                    
    Mn = sum(Mn)+C_c*((h/2)-a/2)
    #return Pn*φ,Mn*φ/12
    return Pn*φ,Mn*φ/12


def find_all_Z_PM(rebar_array,εcu,b,h,fpc,fy=4*ksi):
    PM_points = []
    for Z in np.arange(1,-1000.1,-0.1):
        strains,a,φ = find_strain(Z,rebar_array,εcu)
        forces = find_rebar_force(strains,rebar_array)
        PM_points.append(find_PM(a,b,h,forces,rebar_array,φ))
        
    return PM_points
        
        
        
def plot_truncated_PM(PM_points,rebar_array, b, h, fpc, fy):
    As = 0
    for reb in rebar_array:
        As += reb[0]*rebars[reb[1]]['area']
        
    Max_φPn = 0.8*0.65*(0.85*fpc*((b*h)-(As))+(As*fy))
    
    #PM_points = np.transpose(PM_points)
    #PM_points[0] = np.where(PM_points[0]<Max_φPn,PM_points[0],Max_φPn)
    PM_trunc = np.where(PM_points[0]<Max_φPn,PM_points[0],Max_φPn)
    
    fig, ax = plt.subplots()
    #ax.plot(PM_points[1],PM_points[0])
    ax.plot(PM_points[1],PM_trunc)
    ax.grid()
    ax.set(xlabel='Moment (Kip-ft)', ylabel='Axial Loading (Kip)',
       title='ACI 318 - PM Interaction Diagram')
    plt.show()
    

    
rebar_array = create_rebar_array(num_lay, bar_per_lay, bar_size,CC_to_CL_long)
print(rebar_array)

plot_truncated_PM(np.transpose(find_all_Z_PM(rebar_array,εcu,b,h,fpc,fy=4*ksi)),rebar_array, b, h, fpc, fy)




