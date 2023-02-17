# function to solve a system of differential equations

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

"""
Code to solve a system of differential equations
Taken from: https://scipy-cookbook.readthedocs.io/items/LoktaVolterraTutorial.html

The Lotka-Volterra equations are a pair of first-order, non-linear, differential equations
often used to describe the dynamics of biological systems in which two species interact,
one a predator and the other its prey. The two species are called x and y.
The equations are:

dx/dt = a*x - b*x*y
dy/dt = -c*y + d*b*x*y

where:
a = growth rate of prey when there is no predator
b = death rate of prey due to predation
c = death rate of predator when there is no prey
d = growth rate of predator when there is prey

I have expanded these to include a poaching parameter p
dx/dt = (a - p)*x - b*x*y
dy/dt = -c*y + d*b*x*y

where:
p = poaching rate of outside humans
"""

# Definition of parameters
a = 1.
b = 0.1
c = 1.5
d = 0.75
p = 0.9
def dX_dt(X, t=0):
    #Return the growth rate of fox and rabbit populations
    return [ (a - p)*X[0] -   b*X[0]*X[1] ,
            -(c)*X[1] + d*b*X[0]*X[1] ]

# Define the initial conditions
X0 = [10, 5]

# Define the time range
t = np.linspace(0, 20, 20000)

# Solve the differential equations
X, infodict = odeint(dX_dt, X0, t, full_output=True)

print(infodict['message'])

prey, predator = X.T
f1 = plt.figure()
plt.plot(t, prey, 'r-', label='Prey Population')
plt.plot(t, predator  , 'b-', label='Predator Population')
plt.grid()
plt.xlim(0, 20)
plt.legend(loc='best')
plt.xlabel('Time')
plt.ylabel('Populations')
plt.suptitle('Predator x Pray: Modified Lotka-Volterra Model', fontsize=14)
plt.title("Poaching parameter = {}".format(p), fontsize=10)
f1.savefig('predator_pray.png')