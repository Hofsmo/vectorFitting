import unittest
import numpy as np
import vectorFitting as vf
from numpy import testing


class convolutionTest(unittest.TestCase):
    """Class to perform unit tests on the convolution function"""

    def setUP(self):
        pass

    def test_simple_convolution(self):
        """Test convolution between a unit step function and two poles"""
        t = np.arange(0, 100, 0.01)
        x = (np.ones(len(t))-np.exp(-1000*t))
        waves = vf. windowConv(x, np.array([-1, -2]), t)

        testing.assert_allclose(waves[1],
                                (np.ones(len(t))-np.exp(-2*t))/2, rtol=1e-2, atol=0.01,
                                verbose=True)

if __name__ == '__main__':
    unittest.main()
