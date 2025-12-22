# Functions shared by Lithwick (2012) and Deck & Agol (2016) literatures
import numpy as np
from scipy.integrate import quad
from scipy.special import comb as nCr

# Determine MMR for string input
def get_MMR(MMRstr: str):
    """
    Decomposes string of MMR in form j:j-N into individual j and N as integers.

    For example, '3:1' -> j,N = 3,2 or '5:4' -> j,N = 5,1.

    Moreover, output must also satisfy j > j-N > 0.

    Args:
        MMRstr: String denoting j:j-N mean motion resonance of system e.g. '2:1', '3:2', '3:1', etc. (string)

    Returns:
        j: Number denoting j:j-N mean motion resonance of system (integer)
        N: Number denoting Nth order mean motion resonance of system (integer)
    """
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
    """
    Returns ratio of orbital semi-major axis, computed from Kepler's third law using orbital periods.

    Args:
        innerPeriod: Orbital period of inner planet, in days (float or array)
        outerPeriod: Orbital period of outer planet, in days (float or array)

    Returns:
        alpha: Ratio of orbital semi-major axes. Value is strictly less than 1. (float or array)
    """
    return (innerPeriod / outerPeriod)**(2/3)

# Period calculation from superperiod; Will return two possible solutions
def get_innerPeriods(outerPeriod, superperiod: float, j: int, N: int):
    """
    Returns two possible inner period solutions given outer period and outer TTV superperiod.
    This is computed using the definition of the TTV superperiod (Equation 5 of Lithwick 2012).

    Args:
        outerPeriod: Orbital period of outer planet, in days (float)
        superperiod: TTV superperiod of outer planet, in days (float)

    Returns:
        innerPeriod: Two possible solutions of inner orbital period, in days (array)
    """
    innerPeriod = (j-N) / (j / outerPeriod + 1 / np.array([superperiod, -superperiod]))
    return innerPeriod

def get_outerPeriods(innerPeriod, superperiod: float, j: int, N: int):
    """
    Returns two possible outer period solutions given inner period and inner TTV superperiod.
    This is computed using the definition of the TTV superperiod (Equation 5 of Lithwick 2012).

    Args:
        innerPeriod: Orbital period of inner planet, in days (float)
        superperiod: TTV superperiod of inner planet, in days (float)

    Returns:
        outerPeriod: Two possible solutions of outer orbital period, in days (array)
    """
    outerPeriod = j / ((j-N) / innerPeriod + 1 / np.array([superperiod, -superperiod]))
    return outerPeriod

# Normalized distance to resonance
def get_NormalizedResonanceDistance(innerPeriod, outerPeriod, j: int, N: int):
    """
    Returns "normalized distance to resonance,"
    a unitless parameter denoting how "far" a system is from mean motion resonance
    (Equation 6 of Lithwick 2012).

    Args:
        innerPeriod: Orbital period of inner planet, in days (float or array)
        outerPeriod: Orbital period of outer planet, in days (float or array)
        j: Number denoting j:j-N mean motion resonance of system (integer)
        N: Number denoting Nth order mean motion resonance of system (integer)

    Returns:
        Delta: Normalized distance to resonance (float or array)
    """
    Delta = outerPeriod/innerPeriod * (j-N)/j - 1
    return Delta

# Laplace coefficient and their numerical derivatives
def get_b(alpha, j: int, epsrel=1e-5, method='powerseries'):
    """
    Returns the Laplace coefficient used in three-body perturbation theory of celestial mechanics.
    Since the Laplace coefficient is defined by a nonelementary integral,
    numerical methods must be used to compute it.

    Args:
        alpha: Ratio of orbit semi-major axes, see natsume.common.get_alpha (float or array)
        j: Number denoting j:j-N mean motion resonance of system (integer)
        epsrel: Relative error of computed Laplace coefficient (float, default 1e-5)
        method: 'integrate' or 'powerseries' only (string, default 'powerseries')

        'integrate' computes the Laplace coefficient using numerical integration
        via scipy's scipy.integrate.quad().

        'powerseries' computes the Laplace coefficient using the power series represention as described in
        Appendix C of Marding 2013. Faster than scipy's numerical integration.

    Returns:
        b: The Laplace coefficient (float or array)
    """
    if method == 'integrate': # Numerical integration, ~2x slower
        def func(theta, alpha, j): # Function to integrate
            return np.cos(j*theta) / np.sqrt(1 - 2*alpha*np.cos(theta) + alpha**2)

        if isinstance(alpha, np.ndarray):
            integral = np.array([])
            for a in alpha:
                int, err = quad(func, 0, 2*np.pi, args=(a, j), epsrel=epsrel)
                integral = np.append(integral, int)
        else:
            integral, err = quad(func, 0, 2*np.pi, args=(alpha, j), epsrel=epsrel)

        return integral / np.pi
    
    elif method == 'powerseries': # Power series (Mardling 2013, Appendix C)
        # This method has relative error epsrel = alpha^(2p) where p is the number of iterations
        p_to_epsrel = np.ceil(0.5 * np.log(epsrel) / np.log(alpha))
        p_list = range(round(np.max(p_to_epsrel)) + 1) # p=0 to p=p_to_epsrel.max
        summation = 0
        for p in p_list:
            summation += 2**(1 - 4*p - 2*j) * nCr(2*p, p) * nCr(2*(j+p), j+p) * alpha**(2*p)
        
        laplace_coeff = alpha**j * summation
        return laplace_coeff

    else:
        raise ValueError(f'Laplace coefficient retrieval method {method} does not exist!')

# Central finite difference numerical differentiation (Fornberg 1988)
# All has O(eps^2) error
def get_Db(alpha, j: int, order=1, eps=1e-5):
    """
    Returns the derivatives with respect to alpha of Laplace coefficients using central finite difference
    numerical differentiation as described in Fornberg 1988.

    Args:
        alpha: Ratio of orbit semi-major axes, see natsume.common.get_alpha (float or array)
        j: Number denoting j:j-N mean motion resonance of system (integer)
        order: Order of derivative with respect to alpha. Only supports up to 4th order (integer)
        eps: Epsilon used to compute the Laplace coefficient derivative.
        The numerical differentation yields O(eps^2) error (float, default 1e-5)

    Returns:
        Db: Up to 4th order derivative w.r.t. alpha of the Laplace coefficient (float or array)
    """
    if order == 1:
        return (get_b(alpha + eps, j) - get_b(alpha - eps, j)) / (2 * eps)
    elif order == 2:
        return (get_b(alpha + eps, j) - 2*get_b(alpha, j) + get_b(alpha - eps, j)) / (eps**2)
    elif order == 3:
        return (0.5*get_b(alpha + 2*eps, j) - get_b(alpha + eps, j) + \
                get_b(alpha - eps, j) - 0.5*get_b(alpha - 2*eps, j)) / (eps**3)
    elif order == 4:
        return (get_b(alpha + 2*eps, j) - 4*get_b(alpha + eps, j) + 6*get_b(alpha, j) - \
                4*get_b(alpha - eps, j) + get_b(alpha - 2*eps, j)) / (eps**4)
    else:
        raise ValueError(f'NATSUME does not support numerical differentation of Laplace coefficients at order {order} (N > 4)')
