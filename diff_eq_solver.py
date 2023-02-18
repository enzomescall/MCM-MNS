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
a = 1.0
b = 0.1
c = 1.5
d = 0.75
p = 0.0

X_f0 = np.array([     0. ,  0.])
X_f1 = np.array([ c/(d*b), a/b])

def dX_dt(X, t=0):
    #Return the growth rate of fox and rabbit populations
    return [ (a - p)*X[0] -   b*X[0]*X[1] ,
            -c*X[1] + d*b*X[0]*X[1] ]



# Define the initial conditions
X0 = [10, 5]

# Define the time range
t = np.linspace(0, 20, 20000)

# Solve the differential equations
X, infodict = odeint(dX_dt, X0, t, full_output=True)

prey, predator = X.T

f1 = plt.figure()
plt.plot(t, prey, 'r-', label='Prey 1 Population')
plt.plot(t, predator  , 'b-', label='Predator 1 Population')
plt.grid()
plt.xlim(0, 20)
plt.legend(loc='best')
plt.xlabel('Time')
plt.ylabel('Populations')
plt.suptitle('Predator x Pray: Modified Lotka-Volterra Model', fontsize=14)
plt.title("Poaching parameter = {}".format(p), fontsize=10)
f1.savefig('predator_pray.png')


"""
Now we will plot the trajectories of the system in the (x, y) plane.
We will use the function odeint to solve the system of differential equations.
We will use the function quiver to plot the direction field.
"""

values  = np.linspace(0.3, 0.9, 5)                          # position of X0 between X_f0 and X_f1
vcolors = plt.cm.autumn_r(np.linspace(0.3, 1., len(values)))  # colors for each trajectory

f2 = plt.figure()

# plot trajectories
for v, col in zip(values, vcolors):
    X0 = v * X_f1                     # starting point
    X = odeint( dX_dt, X0, t)         # we don't need infodict here
    plt.plot( X[:,0], X[:,1], lw=3.5*v, color=col, label='X0=(%.f, %.f)' % ( X0[0], X0[1]) )

# define a grid and compute direction at each point
ymax = plt.ylim(ymin=0)[1]                        # get axis limits
xmax = plt.xlim(xmin=0)[1]
nb_points = 20

x = np.linspace(0, xmax, nb_points)
y = np.linspace(0, ymax, nb_points)

X1 , Y1  = np.meshgrid(x, y)                       # create a grid
DX1, DY1 = dX_dt([X1, Y1])                      # compute growth rate on the gridt
M = (np.hypot(DX1, DY1))                           # Norm of the growth rate 
M[ M == 0] = 1.                                 # Avoid zero division errors 
DX1 /= M                                        # Normalize each arrows
DY1 /= M

# Drow direction fields, using matplotlib 's quiver function

plt.title('Trajectories and direction fields')
Q = plt.quiver(X1, Y1, DX1, DY1, M, pivot='mid', cmap=plt.cm.jet)
plt.xlabel('Number of prey')
plt.ylabel('Number of predators')
plt.legend()
plt.grid()
plt.xlim(0, xmax)
plt.ylim(0, ymax)
plt.show()