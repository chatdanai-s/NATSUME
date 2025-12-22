import numpy as np

### Classes
# Complex eccentricities containing the eccentricity and
# longitude of periastron for inner and outer planets (input in deg, kept as rad!) 
class ComplexEccentricities:
    def __init__(self,
                 inner_e=0, inner_periastron=0,
                 outer_e=0, outer_periastron=0):
        """
        Returns complex eccentricities (eccentricity and longitude of periastron) of two planets in the system.

        Used as input in functions requiring eccentricies and longitudes of periastron.

        Args:
            inner_e: Inner orbital eccentricity, must be at least 0 (float)
            inner_periastron: Inner longitude of periastron, in degrees (float)
            outer_e: Outer orbital eccentricity, must be at least 0 (float)
            outer_periastron: Outer longitude of periastron, in degrees (float)

        Returns:
            Complex eccentricities of two planets in the system (class ComplexEccentricities)
        """
        try:
            if (inner_e < 0) or (outer_e < 0):
                raise ValueError(f'inner_e and outer_e in ComplexEccentricities must be greater than or equal to 0.')
            self.inner_e = float(inner_e)
            self.inner_periastron = np.deg2rad(float(inner_periastron))
            self.outer_e = float(outer_e)
            self.outer_periastron = np.deg2rad(float(outer_periastron))

        except (TypeError, ValueError):
            raise TypeError(f'All ComplexEccentricities arguments must be a float or integer.')

    @property
    def arr(self):
        """
        Returns array representation of values in class ComplexEccentricties object.

        Note that the angles are stored and output in rad even if the input is in degrees.

        Args:
            None

        Returns:
            Array object containing inner eccentricity, inner periastron (rad),
            outer eccentricity, and outer periastron (rad), in order
        """
        return np.array([self.inner_e, self.inner_periastron,
                         self.outer_e, self.outer_periastron])

# TTV Sine curve with amplitude in minutes and superperiod in days
class TTVSineCurve:
    def __init__(self, amplitude: float, superperiod: float):
        """
        Returns TTV sine curve representation (Amplitude and Superperiod) of a planet, both in days.
        Does not specify if inner TTV or outer TTV.

        Used as input in functions requiring sinusoidal TTV curves.

        Args:
            amplitude: TTV amplitude in days, must be at least 0 (float)
            superperiod: TTV superperiod in days, must be at least 0 (float)

        Returns:
            TTV sine curve representation of a planet (class TTVSineCurve)
        """
        try:
            if (amplitude <= 0) or (superperiod <= 0):
                raise ValueError(f'Amplitude and superperiod in TTVSineCurve must be greater than 0.')
            self.amplitude = float(amplitude)
            self.superperiod = float(superperiod)

        except (TypeError, ValueError):
            raise TypeError(f'All TTVSineCurve arguments must be a float or integer.')
        
    @property
    def arr(self):
        """
        Returns array representation of values in class TTVSineCurve object.

        Args:
            None

        Returns:
            Array object containing TTV amplitude and TTV superperiod, in order
        """
        return np.array([self.amplitude, self.superperiod])
