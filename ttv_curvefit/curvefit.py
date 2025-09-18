# Curve fits TTVs and returns TTV amplitude and superperiod
# Make sure the initial period guess is close to actual period! Else the fitting doesn't work
# Lomb-Scargle periodogram not used because superperiod can already be estimated from TTVFast inputs
import numpy as np
from scipy.optimize import curve_fit

def return_ttv(epochs, midtransits):
    # TTVFast output
    n = np.array(epochs)
    tn = np.array(midtransits)

    # Constant-period model
    def linear_model(n, P, t0):
        return t0 + n * P

    # Finds average period
    p0 = [tn[-2] - tn[-1],
         0.5 * (tn[0] + tn[-1])]  # Initial guess for P and t0
    popt_lin, pcov_lin = curve_fit(f=linear_model, xdata=n, ydata=tn, p0=p0)
    perr_lin = np.sqrt(np.diag(pcov_lin))
    
    # Extract variables
    P, t0 = popt_lin
    Perr, t0err = perr_lin

    # Get TTV signal (O minus C per time)
    time = tn
    ttv = tn - linear_model(n, P, t0)
    ttv_err = np.sqrt(t0err**2 + (n*Perr)**2)

    return time, ttv, ttv_err

def ttv_model(t, A, B, V, Pttv, phase):
    return A + B*t + V * np.sin(2*np.pi/Pttv * t + phase)
    
def characterize_ttv(time, ttv, ttv_err, expected_period):
    # Fit sine curve to TTV
    V_guess = (ttv.max() - ttv.min()) / 2
    Pttv_guess = expected_period
    phase_guess = (np.pi/2 - 2*np.pi/Pttv_guess * time[np.argmax(ttv)]) % (2 * np.pi)

    p0 = [0, 0, V_guess, Pttv_guess, phase_guess]  # Initial guess for parameters
    popt_ttv, pcov_ttv = curve_fit(f=ttv_model, xdata=time, ydata=ttv, p0=p0,
                                   bounds=([-np.inf, -np.inf, 0, 0, -np.pi],
                                           [np.inf, np.inf, 2*V_guess, np.inf, 3*np.pi]))
    perr_ttv = np.sqrt(np.diag(pcov_ttv))

    # Acquire reduced chi-squares
    residuals = ttv - ttv_model(time, *popt_ttv)
    red_chi2 = np.sum((residuals/ttv_err)**2) / (len(ttv) - len(popt_ttv))

    # Return [A, B, V, Pttv, phase], their errors, and red-chi2 score
    return popt_ttv, perr_ttv, red_chi2
