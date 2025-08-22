# Functions shared by Lithwick (2012) and Deck & Agol (2016) literatures
import numpy as np
from scipy.integrate import quad

# Determine MMR for string input
def get_MMR(MMRstr: str):
    try: 
        j, i = MMRstr.split(':')
        jj, ii = int(j), int(i)
        if (jj == float(j)) and (ii == float(i)) and (jj > ii) and (ii > 0):
            N = jj - ii
            return jj, N
        else:
            raise ValueError(f"MMR argument must be in 'j:i' format, where i, j are integers such that j > i > 0.")

    except ValueError:
        raise ValueError(f"MMR argument must be in 'j:i' format, where i, j are integers such that j > i > 0.")
    
# Semi-major axis ratio
def get_alpha(innerPeriod, outerPeriod):
    return (innerPeriod / outerPeriod)**(2/3)

# Period calculation from superperiod; Will return two possible solutions
def get_innerPeriods(outerPeriod, superperiod: float, j: int, N: int):
    innerPeriod = (j-N) / (j / outerPeriod + 1 / np.array([superperiod, -superperiod]))
    return innerPeriod

def get_outerPeriods(innerPeriod, superperiod: float, j: int, N: int):
    outerPeriod = j / ((j-N) / innerPeriod + 1 / np.array([superperiod, -superperiod]))
    return outerPeriod

def get_NormalizedResonanceDistance(innerPeriod, outerPeriod, j: int, N: int):
    Delta = outerPeriod/innerPeriod * (j-N)/j - 1
    return Delta

# Laplace coefficient and their numerical derivatives
def get_b(alpha, j: int):
    def func(theta, alpha, j): # Function to integrate
        return np.cos(j*theta) / np.sqrt(1 - 2*alpha*np.cos(theta) + alpha**2)

    if isinstance(alpha, np.ndarray):
        integral = np.array([])
        for a in alpha:
            int, err = quad(func, 0, 2*np.pi, args=(a, j))
            integral = np.append(integral, int)
    else:
        integral, err = quad(func, 0, 2*np.pi, args=(alpha, j))
    
    return integral / np.pi

def get_db_da(alpha, j: int, eps=1e-6): # O(eps^2) error
    return (get_b(alpha + eps, j) - get_b(alpha - eps, j)) / (2 * eps)

def get_d2b_da2(alpha, j: int, eps=1e-6): # O(eps^2) error
    return (get_b(alpha + eps, j) - 2*get_b(alpha, j) + get_b(alpha - eps, j)) / (eps**2)
