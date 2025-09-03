# Curve fits TTVs and returns TTV amplitude and superperiod
import numpy as np
from scipy.optimize import curve_fit

def characterize_ttv(epochs, midtransits, expected_period, return_graph=False):
    # TTVFast output
    n = np.array(epochs)
    tn = np.array(midtransits)

    # Sinusoidal TTV model from mid-transit times
    def linear_model(n, P, t0):
        return t0 + n * P

    def ttv_model(t, V, Pttv, phase):
        return V * np.sin(2*np.pi/Pttv * t + phase)

    # Finds average period
    p0 = [tn[1] - tn[0],
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

    # Fit sine curve to TTV
    V_guess = (ttv.max() - ttv.min()) / 2
    Pttv_guess = expected_period
    phase_guess = (np.pi/2 - 2*np.pi/Pttv_guess * time[np.argmax(ttv)]) % (2 * np.pi)

    p0 = [V_guess, Pttv_guess, phase_guess]  # Initial guess for parameters
    popt_ttv, pcov_ttv = curve_fit(f=ttv_model, xdata=time, ydata=ttv, p0=p0,
                                   bounds=([0, 0, -np.pi],[2*V_guess, np.inf, 3*np.pi]))
    perr_ttv = np.sqrt(np.diag(pcov_ttv))

    # Acquire reduced chi-squares
    residuals = ttv - ttv_model(time, *popt_ttv)
    red_chi2 = np.sum((residuals/ttv_err)**2) / (len(ttv) - len(popt_ttv))

    # Append parameters 
    param_final = np.append(popt_lin, popt_ttv)
    param_f_err = np.append(perr_lin, perr_ttv)

    # TBA: Graph to check fit by eye
    if return_graph == True:
        # TTV part
        # Residue part
        # Labels
        pass

    # Return [t0, P, V, Pttv, phase], their errors, and red-chi2 score
    return param_final, param_f_err, red_chi2
