# This file contains the function HOC, used to randomly generate the 3d coordinates for Hafnium Carbide.
import numpy as np
from random import random

'''
FIX COMMENTS
Function to define the probability of sampling carbon at layer z given a linear change between the outer layers of the 
Hafnium carbide material, and the middle layer(s), given the expected percent of Carbon versus Oxygen overall and at the outermost layers. 
The assumption is that the outer layers have more Oxygen or Carbon, and the inner layers have a greater mixture of the 
two atoms.
The change from outer edge layer to inner mixed layers is linear, with the rate of change determined by the number of 
layers in the boundary.
Input z, integer, the specific layer of the material being considered with index starting at 0. This 
         function can be used to define a lambda functional with z as the variable.
Input maxz, even integer, the number of layers in the material, should be twice the number of Hafnium in each vertical column.
Input carbonPercent, float in [0,1], the expected average percent of overall carbon in the material out of the total Carbon or Oxygen atoms.
Input outerPercent,  list of floats in [0,1], the average percent of carbon in the topmost outer layers is outerPercent[0], and 
         the average percent of carbon at the bottom outer layers is outerPercent[1]. The default is outerPercent  = [0, 1].
Input outerLayer, list of two integers. outerLayer[0] is the number of layers at the top with percent carbon = outerPercent[0],
         and outerLayer[1] is the number of layers at the bottom of the material with percent carbon = outerPercent[1].
         The default is [1,1].
Input boundary, list of two integers. boundary[0] is the number of layers below outerLayer[0] that the percent of Carbon 
         linearly transitions from outerPercent[0] to carbonPercent over. boundary[1] is likewise the number of layers above 
         outerLayer[1] that it takes to linearly transition between. The default is that the two layers exactly in the middle
         of the material return carbonPercent, and all other unaccounted for layers linearly transition between outer and the middle.
Output, the expected percent of Carbon at layer z (given the inputs and linear transition assumption).  
'''
def zvariation(z, maxz, carbonPercent, boundaryPercent = [1,0], boundary = [1,1], transition = 0):
    if transition == 0:
        transition[0] = maxz/2-1
        transition[1] = transition[0] - boundary[0]
        transition[0] = transition[0] - boundary[1]
    else:
        if isinstance(transition[0], (int, float)) & isinstance(transition[1], (int, float)) & isinstance(boundary[0], (int, float)) & isinstance(boundary[1], (int, float)):
            if sum(transition)+sum(boundary) <= maxz:
                pass
            else:
                print("ERROR: The number of boundary and transition layers must be less than or equal to maxz.")
        else:
            print("ERROR: boundary and transition must each be an array with two integers")
    if z < boundary[0]:
        return boundaryPercent[0]
    elif z < boundary[0]+transition[0]:
        return boundaryPercent[0] - (boundaryPercent[0]-carbonPercent)/transition[0]*z
    elif maxz - z <= boundary[1]:
        if maxz - z < 1:
            print("WARNING: The vertical distribution you are using assumes fewer total layers.")
        return boundaryPercent[1]
    elif maxz - z <= boundary[1] + transition[1]:
        return carbonPercent - (boundaryPercent[1]-carbonPercent)/transition[1]*z
    else:
        return carbonPercent

'''
 hfoc is used to randomly generate the coordinates for Oxygen and Carbon molecules within a Hafnium Carbide Molecule
 input dims, array of 3 integers, number of atoms in the x direction, y direction, then z direction
 input verticalDist, function with z as input, describes how the probability of finding a Carbon atom changes for each layer of z.
                     The default value is to use a uniform distribution of O and C for all layers.
                     See the function zvariation in this same module for more options.

 output HfOCcoordinates, numpy array of floats, with coordinates of Hafnium, followed by coordinates of Oxygen, followed 
                     coordinates of Carbon atoms in the randomly generated material.
                     It has 4*dim[1]*dim[2]*dim[3] atoms, each represented by a row, and three columns for the [x, y, z]
                     coordinates respectively. The bottom left of the material is at [0 0 0], and surrounding atoms are added 
                     in the positive coordinate directions. 
 output AtomNumbers, integer array, the number of Hafnium atoms in the first rows of HOCcoordinates,
                     the number of Oxygen atoms in the middle rows of HOCcoordinates, and
                     the number of Carbon atoms in the last rows of HOCcoordinates.
'''
def hfoc(dims, verticalDist = 0):
    #Parse inputs
    if len(dims) != 3:
        print("the first input must be an array with three positive integers")
        return
    elif dims[0] < 1:
        dims[0] = 1
        print("Warning: x dimension set to 1")
    elif dims[1] < 1:
        dims[1] = 1
        print("Warning: y dimension set to 1")
    elif dims[2] < 1:
        dims[2] = 1
        print("Warning: z dimension set to 1")
    #
    if verticalDist == 0:
        verticalDist = lambda z: 0.5
    #
    #Define place to save outputs
    hf = np.array([[-1, -1, -1]])
    o = np.array([[-1, -1, -1]])
    c = np.array([[-1, -1, -1]])
    #
    #Begin generating the molecule
    for z in np.arange(0, dims[2]):
        cutOff0 = verticalDist(2*z)
        cutOff1 = verticalDist(2*z+1)
        for y in np.arange(0, dims[1]):
            for x in np.arange(0, dims[0]):
                hf = np.concatenate((hf, np.array([[2*x,2*y,2*z]]), np.array([[2*x+1, 2*y+1, 2*z]]), np.array([[2*x+1, 2*y, 2*z+1]]), np.array([[2*x, 2*y+1, 2*z+1]])), axis=0)
                #randomly generate 4 numbers in zero 1
                m = [random(), random(), random(), random()]
                #save the next coordinates as O or C according to the m
                for ii in [0, 1]:
                    if m[ii] < cutOff0:
                        c = np.concatenate((c,np.array([[2*x+1-ii,2*y+ii,2*z]])), axis = 0)
                    else:
                        o = np.concatenate((o,np.array([[2*x+1-ii, 2*y+ii, 2*z]])), axis=0)
                    if m[ii+2] < cutOff1:
                        c = np.concatenate((c, np.array([[2*x+ii,2*y+ii,2*z+1]])), axis=0)
                    else:
                        o = np.concatenate((o, np.array([[2*x+ii, 2*y+ii, 2*z+1]])), axis=0)
    # Remove the first filler rows
    hf = hf[1:, :]
    c = c[1:, :]
    o = o[1:, :]
    HfOCcoordinates = 2.305*np.concatenate((hf, o, c), axis=0)
    AtomNumbers = [len(hf), len(o), len(c)]
    return HfOCcoordinates, AtomNumbers



