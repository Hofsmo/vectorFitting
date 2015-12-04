import numpy as np
from findPoles import findPoles

<<<<<<< HEAD
t = np.arange(0,20,0.02)

x = 2*np.ones((1,len(t)))

y = np.exp(-t)

p = np.array([0.8, -3])

poles =  findPoles(x,y,p,t)
=======
def vectorFitting(x, y, t, tol=1e-5, iMax)
    """Function that performs time domain vector fitting

    Input:
        x is the input to the system
        y is the system's response to the input x
        initPoles is a numpy column vector containing the inital poles
        t is the time vector"""
    for i in np.arange(iMax):
        if np.linalg.norm(poles-initPoles) < tol or iterator >= 100:
            poles =  findPoles(x,y,p,t)
