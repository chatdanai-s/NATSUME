# Curve fits TTVs and returns TTV amplitude and superperiod
import numpy as np
from scipy.optimize import curve_fit

def characterize_ttv(epochs, midtransits, expected_cycles=1):
    # TTVFast output
    n = np.array(epochs)
    tn = np.array(midtransits)

    # Sinusoidal TTV model from mid-transit times
    def linear_model(n, P, t0):
        return t0 + n * P

    def ttv_model(t, A, B, V, Pttv, phase):
        return A + B*t + V * np.sin(2*np.pi/Pttv * t + phase)

    # Finds average period
    p0 = [tn[1] - tn[0],
         0.5 * (tn[0] + tn[-1])]  # Initial guess for P and t0
    popt_lin, pcov_lin = curve_fit(f=linear_model, xdata=n, ydata=tn, p0=p0)
    perr_lin = np.sqrt(np.diag(pcov_lin))

    P, t0 = popt_lin
    Perr, t0err = perr_lin

    # Get TTV signal (O minus C per time)
    time = n * P
    time_err = np.abs(n) * Perr
    ttv = tn - linear_model(n, *popt_lin)
    ttv_err = None # Not yet implemented

    # Fit sine curve to TTV
    p0 = [0, 0,
          (ttv.max() - ttv.min()) / 2,
          (time.max() - time.min()) / expected_cycles,
          0] # Initial guess for parameters
    popt_ttv, pcov_ttv = curve_fit(f=ttv_model, xdata=time, ydata=ttv, p0=p0) # sigma arg soonTM
    perr_ttv = np.sqrt(np.diag(pcov_ttv))

    # Return A, B, V, Pttv, phase and their errors
    return popt_ttv, perr_ttv
