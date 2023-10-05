import matplotlib.pyplot as plt
import numpy as np
import h5py
from load_and_taylor import load_data_from_mat, fit_polynomial, create_taylor_approximation, plot_approximation, k_lambda
from scipy.optimize import fsolve

def lsli(lamp, NL, as_, coeff):
    # Define the function for which we seek the root
    def DK(lamp, lams, coeff):
        LAMI = 1. / (2. / lamp - 1. / lams)
        
        kp = k_lambda(lamp, coeff)
        ks = k_lambda(lams, coeff)
        ki = k_lambda(LAMI, coeff)
        
        return kp + kp - ks - ki - 1e-6


    # Use fsolve to find the root
    lams = fsolve(DK, as_, args=(lamp, coeff, NL), xtol=1e-20, maxfev=10000)[0]
    lami = (2 * 3.14159 * 3e14) / ((2 * 3.14159 * 3e14) / lamp + (2 * 3.14159 * 3e14) / lamp - (2 * 3.14159 * 3e14) / lams)

    return lams, lami



path = '/home/jay/repos/F3001C_Reto/CodigoFinal/Modos/Modes/SweepOverlapTE/Waveguide1000_550_1580.mat'

neff_values, lambda_values = load_data_from_mat(path)

# Fit polynomial
coefficients = fit_polynomial(lambda_values, neff_values)

# Create Taylor approximation
polynomial = np.poly1d(coefficients)
center_lambda = (lambda_values.min() + lambda_values.max()) / 2
taylor_polynomial = create_taylor_approximation(polynomial, center_lambda)

# Evaluate polynomials
lambda_range = np.linspace(lambda_values.min(), lambda_values.max(), 500)
original_values = np.polyval(coefficients, lambda_range)
taylor_values = taylor_polynomial(lambda_range - center_lambda)

# Plot Neff approximation
plot_approximation(lambda_values, neff_values, lambda_range, original_values, taylor_values, 'Neff', 'Taylor Polynomial Approximation of Neff')

L = 0.1e6
sigma = 0.1e12
lamp0 = 0.762
NL = 0

omgp0 = (2*np.pi*3e14)/lamp0
omgp = np.linspace(omgp0-3*sigma, omgp0+3*sigma,100)

a_s = 0.6

[lams0, lami0] = lsli(lamp0, NL, a_s, coefficients)

omgs0 = (2*np.pi*3e14)/lams0
omgi0 = (2*np.pi*3e14)/lami0

dw = 24e12
Ns = 100



oms = np.linspace(omgs0 - dw, omgs0 + dw, Ns)


