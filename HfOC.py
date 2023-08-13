# This file contains the function HOC, used to randomly generate the 3d coordinates for Hafnium Carbide.
import numpy as np

'''
Function to define the probability of sampling carbon at layer z given a linear change between the vertical layers of the 
Hafnia Carbide material. The assumption is that the outer layers have more Oxygen or Carbon, and the inner layers have a greater mixture of the 
two atoms. The change from outer edge layer to inner mixed layers is linear, where the slopes depend on the function inputs.
Input z, integer, the specific layer of the material being considered with index starting at 0. This function is designed 
         for defining a lambda functional with z as the variable.
Input maxz, even integer, the number of layers in the material, equal to twice the number of Hafnium in each vertical column.
Input carbonPercent, float in [0,1], the expected average percent of carbon versus oxygen atoms in the central layers of the material.
Input boundaryPercent,  list of 2 floats in [0,1], the average percent of carbon in the bottom boundary layers is outerPercent[0], and 
         the average percent of carbon at the top boundary layers is outerPercent[1]. The default is boundaryPercent  = [1, 0].
Input boundary, list of two integers. boundary[0] is the number of layers in the bottom boundary layers,
         and boundary[1] is the number of layers in the top boundary layers of the material. The default is boundary  = [1,1].
Input transition, list of two integers. transition[0] is the number of layers above the bottom boundary layers that the percent of Carbon 
         linearly transitions from boundaryPercent[0] to carbonPercent over. transition[1] is likewise the number of layers below
         boundary[1] that the probability of choosing Carbon linearly transitions between the two percents. 
         The default is that the two layers exactly in the middle of the material return carbonPercent, and all other 
         unaccounted for layers linearly transition between boundary and those middle layers.
Output, the expected percent of Carbon at layer z (given the inputs and linear transition assumption). 
*Note that the overall percent of Carbon in the material is not equal to carbonPercent, but to the sum of outputs for all z 
 divided by maxz
Author: Natalie Wellen 
'''
def zvariation(z, maxz, carbonPercent, boundaryPercent = [1,0], boundary = [1,1], transition = [-1,-1]):
    if transition[0] == -1:
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
        return boundaryPercent[0] - (boundaryPercent[0]-carbonPercent)/transition[0]*(z+1-boundary[0])
    elif maxz - z <= boundary[1]:
        if maxz - z < 1:
            print("WARNING: The vertical distribution you are using assumes fewer total layers.")
        return boundaryPercent[1]
    elif maxz - z <= boundary[1] + transition[1]:
        return boundaryPercent[1] + (carbonPercent-boundaryPercent[1])/transition[1]*(maxz-z-boundary[1])
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
Author: Natalie Wellen
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
                m = np.random.rand(2,2)
                #save the next coordinates as O or C according to the m
                for ii in [0, 1]:
                    if m[0, ii] < cutOff0:
                        c = np.concatenate((c,np.array([[2*x+1-ii,2*y+ii,2*z]])), axis = 0)
                    else:
                        o = np.concatenate((o,np.array([[2*x+1-ii, 2*y+ii, 2*z]])), axis=0)
                    if m[1, ii] < cutOff1:
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

#Future work: use the multiprocessing module to parallelize the for loops.

