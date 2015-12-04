import numpy as np
import vectorFitting as vf

t = np.arange(0, 100, 0.01)
x = np.ones((len(t), 1))

waves = vf. windowConv(x, t, np.array([-1, -2]))
