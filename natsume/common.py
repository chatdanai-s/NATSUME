# Functions shared by Lithwick (2012) and Deck & Agol (2016) literatures
import numpy as np
from scipy.integrate import quad


# Determine MMR for string input
def MMR(MMRstr: str):
    try: 
        j, i = MMRstr.split(':')
        jj, ii = int(j), int(i)
        if (jj == float(j)) and (ii == float(i)) and (jj > ii) and (ii > 0):
            N = jj - ii
            return jj, N
        else:
            raise ValueError(f"MMR argument must be in 'j:i' format, where i and j are positive integers such that j > i.")

    except ValueError:
        raise ValueError(f"MMR argument must be in 'j:i' format, where i and j are positive integers such that j > i.")
    
# Semi-major axis ratio
def alpha(innerPeriod, outerPeriod):
    return (innerPeriod / outerPeriod)**(2/3)

# Period calculation from superperiod; Will return two possible solutions
def InnerPeriods(outerPeriod, superperiod, j: int, N: int):
    innerPeriod = (j-N) / (j / outerPeriod + 1 / np.array([superperiod, -superperiod]))
    return innerPeriod

def OuterPeriods(innerPeriod, superperiod: float, j: int, N: int):
    outerPeriod = j / ((j-N) / innerPeriod + 1 / np.array([superperiod, -superperiod]))
    return outerPeriod

def NormalizedResonanceDistance(innerPeriod, outerPeriod, j: int, N: int):
    Delta = outerPeriod/innerPeriod * (j-N)/j - 1
    return Delta

# Laplace coefficient and their numerical derivatives
def b(alpha, j: int):
    def func(theta, alpha, j):
        return np.cos(j*theta) / np.sqrt(1 - 2*alpha*np.cos(theta) + alpha**2)

    integral, err = quad(func, 0, 2*np.pi, args=(alpha, j))
    return integral / np.pi

def db_da(alpha, j: int, eps=1e-6): # O(eps^2) error
    return (b(alpha + eps, j) - b(alpha - eps, j)) / (2 * eps)

def d2b_da2(alpha, j: int, eps=1e-6): # O(eps^2) error
    return (b(alpha + eps, j) - 2*b(alpha, j) + b(alpha - eps, j)) / (eps**2)
