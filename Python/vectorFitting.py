import numpy as np


class fitVector:
    """Class used for performing vetor fitting"""
    def __init__(self, x, y, t, tol=1e-5, iMax=100):
        """Constructor for the fitVector class

        Input:
            x: is the input to the system
            y: is the system's response to the input x
            t: is the time vector
            tol: is the decired precision of the fitting
            iMax: is the maximum number of allowed iterations"""

        self.x = x
        self.y = y
        self.t = t
        self.tol = tol
        self.iMax = iMax

    def fitVector(self, initPoles):
        """Function that performs time domain vector fitting

        Input:
            initPoles is a numpy column vector containing the inital poles

        Output:

                """
        poles = 2*initPoles

        for i in np.arange(self.iMax):
            if np.linalg.norm(poles-initPoles) < self.tol:
                poles = initPoles
                initPoles = self.findPoles(poles)
            else:
                return initPoles

        return []

    def findPoles(self, initPoles):
        """Function that performs vector fitting

        Input:
            initPoles is a numpy column vector containing the initial poles
        Output:
            poles the poles in the system
            """

        # Initialize variables
        n = len(initPoles)  # The number of poles in the system
        timesteps = len(self.t)  # The number of timesteps in the input data

        xn = np.empty((n, timesteps), np.complex)
        yn = np.empty((n, timesteps), np.complex)

        # Create waveforms from convolution
        xn = windowConv(self.x, initPoles, self.t)
        yn = windowConv(self.y, initPoles, self.t)

        # Prepare vectors for least squares
        A = np.concatenate((np.reshape(self.x, (timesteps, 1)),
                            np.ones((timesteps, 1))*2,
                            np.reshape(xn, (timesteps, n)),
                            np.reshape(-yn, (timesteps, n))), 1)
        tempy = np.reshape(self.y, (timesteps, 1))

        # Solve the system using least squares
        sol = np.linalg.lstsq(A, tempy, rcond=1e-4)

        # Find kn in the solution vector sol
        kn = sol[0][len(sol[0])-n:len(sol[0])]
        poles = findZeros(kn, initPoles)

        # Flip unstable poles
        poles.real[poles.real < 0] = -poles.real[poles.real < 0]


def findZeros(kn, qn):
    """Function that finds the zeros in the fit function sigma"""

    return np.linalg.eig(np.diag(qn)-np.ones((len(kn), 1))*kn)


def windowConv(x, poles, t):
    """Calculate the convolution for each poles

    Input:
        x: Is the data to convolve with the poles
        t: Is the time vector
        poles: Are the poles

    Output:
        waves: The waves resulting from the convolution"""

    timesteps = len(t)

    waves = np.array(list(map((lambda pole: np.convolve(
        np.exp(pole*t), x)), poles)))

    return waves[:, 0:timesteps]*(t[2]-t[1])
