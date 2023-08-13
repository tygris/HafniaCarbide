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
if test_hf.all() != 0:
    print("Test 1 failed, Hafnium are in the wrong location.")
test_o = a[b[0]:(b[0]+b[1]), :] - oxygen_locations
if test_o.all() != 0:
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
if test_hf.all() != 0:
    print("Test 2 failed, Hafnium are in the wrong location.")

#Test 3, check that the average values are approximately correct for 50% Carbon versus Oxygen
carbon = 0
enum = 100
for ii in np.arange(enum):
    a, b = hf.hfoc([10,10,10])
    carbon += b[2]
    if b[1] + b[2] != b[0]:
        print("Test failed, the number of Hafnium does not equal the number of Carbon and Oxygen combined.")
    if sum(b) != 20**3:
        print("Test failed, the number of atoms in the material is not correct.")
cp = carbon/enum/10**3/4
if cp > 0.6 or cp < 0.4:
    print("Test 3 failed, the percent of Carbon varies a lot from the input 50%.")

#The next set of tests are to test the helper function zvariation
#The tests are done visually by plotting the results

#Test 4, check the linear function for z
zvar = lambda z: hf.zvariation(z, 14, .5, [1,0], [2,2], [3,3])
x = np.arange(14)
y = -1*np.ones(14)
for ii in x:
    y[ii] = zvar(ii)
plt.plot(x,y)
plt.title("Test 4")

#Test 5, test that the default value of transition works
zvar = lambda z: hf.zvariation(z, 100, .5, [1,0], [10,10])
x = np.arange(100)
y = -1*np.ones(100)
for ii in x:
    y[ii] = zvar(ii)
plt.figure()
plt.plot(x,y)
plt.title("Test 5")

#Test 6, test for different overall carbon percentages
zvar = lambda z: hf.zvariation(z, 100, .4, [1,0], [10,10], [25, 25])
x = np.arange(100)
y = -1*np.ones(100)
for ii in x:
    y[ii] = zvar(ii)
plt.figure()
plt.plot(x,y)
plt.title("Test 6")

#Test 7, test for different transition and boundary lengths on each side
zvar = lambda z: hf.zvariation(z, 100, .7, [1,0], [20,10], [10, 40])
x = np.arange(100)
y = -1*np.ones(100)
for ii in x:
    y[ii] = zvar(ii)
plt.figure()
plt.plot(x,y)
plt.title("Test 7")