import HfOC as hf
import numpy as np
import matplotlib.pyplot as plt

#Test 1, all oxygen and 4 hafnium
zvar = lambda z: 0
dims = [1, 1, 1]
hafnium_locations = np.array([[0., 0., 0.],
                              [2.305, 2.305, 0.],
                              [2.305, 0., 2.305],
                              [0., 2.305, 2.305]])
oxygen_locations = np.array([[2.305, 0., 0.],
                             [0., 0., 2.305],
                             [0., 2.305, 0.],
                             [2.305, 2.305, 2.305]])
a,b = hf.hfoc(dims, zvar)
if b[0]!=4:
    print("Test 1 failed, there is the wrong amount of Hafnium.")
if b[1]!=4:
    print("Test 1 failed, there is the wrong amount of Oxygen.")
test_hf = a[0:b[0], :] - hafnium_locations
if test_hf.all() != 0
    print("Test 1 failed, Hafnium are in the wrong location.")
test_o = a[b[0]:(b[0]+b[1]), :] - oxygen_locations
if test_o.all() != 0
    print("Test 1 failed, Oxygen are in the wrong location.")

#test 2, all oxygen and 32 Hafnium
hafnium_locations = np.array([[0.   , 0.   , 0.   ],
       [2.305, 2.305, 0.   ],
       [2.305, 0.   , 2.305],
       [0.   , 2.305, 2.305],
       [4.61 , 0.   , 0.   ],
       [6.915, 2.305, 0.   ],
       [6.915, 0.   , 2.305],
       [4.61 , 2.305, 2.305],
       [0.   , 4.61 , 0.   ],
       [2.305, 6.915, 0.   ],
       [2.305, 4.61 , 2.305],
       [0.   , 6.915, 2.305],
       [4.61 , 4.61 , 0.   ],
       [6.915, 6.915, 0.   ],
       [6.915, 4.61 , 2.305],
       [4.61 , 6.915, 2.305],
       [0.   , 0.   , 4.61 ],
       [2.305, 2.305, 4.61 ],
       [2.305, 0.   , 6.915],
       [0.   , 2.305, 6.915],
       [4.61 , 0.   , 4.61 ],
       [6.915, 2.305, 4.61 ],
       [6.915, 0.   , 6.915],
       [4.61 , 2.305, 6.915],
       [0.   , 4.61 , 4.61 ],
       [2.305, 6.915, 4.61 ],
       [2.305, 4.61 , 6.915],
       [0.   , 6.915, 6.915],
       [4.61 , 4.61 , 4.61 ],
       [6.915, 6.915, 4.61 ],
       [6.915, 4.61 , 6.915],
       [4.61 , 6.915, 6.915]])
a, b = hf.hfoc([2,2,2], zvar)
if b[0]!=32:
    print("Test 2 failed, there is the wrong amount of Hafnium.")
if b[1]!=32:
    print("Test 2 failed, there is the wrong amount of Oxygen.")
test_hf = a[0:b[0], :] - hafnium_locations
if test_hf.all() != 0
    print("Test 2 failed, Hafnium are in the wrong location.")

#Test 3, check that the average values are approximately correct for 50% Carbon versus Oxygen
carbon = 0
enum = 100
for ii in np.arange(enum):
    a, b = HfOC.hfoc([10,10,10])
    carbon += b[2]
    if b[1] + b[2] != b[0]:
        print("Test failed, the number of Hafnium does not equal the number of Carbon and Oxygen combined.")
    if sum(b) != 20**3:
        print("Test failed, the number of atoms in the material is not correct.")
cp = carbon/enum/10**3/4
if cp > 0.6 or cp < 0.4:
    print("Test 3 failed, the percent of Carbon varies a lot from the input 50%.")


#Test 4, check the linear function for z
zvar = lambda z: zvariation(z, 10, .5, boundaryPercent = [1,0], boundary = [2,2], transition = [3,3])


