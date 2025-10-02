# Curve fits TTVs and returns TTV amplitude and superperiod
# Make sure the initial period guess is close to actual period! Else the fitting doesn't work
# Lomb-Scargle periodogram not used because superperiod can already be estimated from TTVFast inputs
import numpy as np
from scipy.optimize import curve_fit
import lmfit

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
    
def characterize_ttv(time, ttv, ttv_err, expected_period, method='scipy'):
    # Fit sine curve to TTV
    V_guess = (ttv.max() - ttv.min()) / 2
    Pttv_guess = expected_period
    phase_guess = (np.pi/2 - 2*np.pi/Pttv_guess * time[np.argmax(ttv)]) % (2 * np.pi)

    p0 = [0, 0, V_guess, Pttv_guess, phase_guess]  # Initial guess for parameters

    if method == 'scipy':
        popt_ttv, pcov_ttv = curve_fit(f=ttv_model, xdata=time, ydata=ttv, p0=p0,
                                    bounds=([-np.inf, -np.inf, 0, 0, -np.pi],
                                            [np.inf, np.inf, 2*V_guess, np.inf, 3*np.pi]))
        perr_ttv = np.sqrt(np.diag(pcov_ttv))

        # Predictions
        ttv_pred = ttv_model(time, *popt_ttv)

        # R² manually
        ss_res = np.sum((ttv - ttv_pred)**2)
        ss_tot = np.sum((ttv - np.mean(ttv))**2)
        r2 = 1 - (ss_res / ss_tot)

    if method == 'lmfit':
        model = lmfit.Model(ttv_model)
        
        # Fitting part
        result = model.fit(ttv, t=time, A=0, B=0, V=V_guess, Pttv=Pttv_guess, phase=phase_guess)

        # Extract values and uncertainties
        popt_ttv = [result.params[p].value for p in ['A', 'B', 'V', 'Pttv', 'phase']]
        perr_ttv = [result.params[p].stderr for p in ['A', 'B', 'V', 'Pttv', 'phase']]

        # Reduced chi-square
        r2 = result.rsquared

    # Return [A, B, V, Pttv, phase], their errors, and red-chi2 score
    return popt_ttv, perr_ttv, r2
