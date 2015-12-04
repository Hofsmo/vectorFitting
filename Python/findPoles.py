import numpy as np
import scipy.signal as sig


def findPoles(x, y, initPoles, t, tol=1e-3, iterator=0):
    """Function that performs vector fitting

    Input:
        x is the input data to the system
        y is the response to the input x
        initPoles is a numpy column vector containing the initial poles
        t is the time vector.
    Output:
        poles the poles in the system
        """

    # Initialize variables
    n = len(initPoles)  # The number of poles in the system
    timesteps = len(t)  # The number of timesteps in the input data

    xn = np.empty((n, timesteps), np.complex)
    yn = np.empty((n, timesteps), np.complex)

    # The vectors used in the convolution should have the same dimensions
    if x.shape != t.shape:
        x = np.reshape(x, t.shape)

    if y.shape != t.shape:
        y = np.reshape(y, t.shape)

    dt = t[1] - t[0]

    # Create waveforms from convolution
    for i in np.arange(n):
        xn[i] = np.convolve(np.exp(initPoles[i]*t), x)[0:timesteps]*dt
        yn[i] = np.convolve(np.exp(initPoles[i]*t), y)[0:timesteps]*dt

    # Prepare vectors for least squares
    A = np.concatenate((np.reshape(x, (timesteps, 1)),
                        np.ones((timesteps, 1))*2,
                        np.reshape(xn, (timesteps, n)),
                        np.reshape(-yn, (timesteps, n))), 1)
    y = np.reshape(y, (timesteps, 1))

    # Solve the system using least squares
    sol = np.linalg.lstsq(A, y, rcond=1e-4)

    # Find kn in the solution vector sol
    kn = sol[0][len(sol[0])-n:len(sol[0])]
    poles = findZeros(kn, initPoles)



def findZeros(kn, qn):
    """Function that finds the zeros in the fit function sigma"""

    [a, b] = sig.invres(kn, qn, np.array([1]))

    return np.roots(a)
