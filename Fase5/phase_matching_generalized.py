import matplotlib.pyplot as plt
import numpy as np
import h5py
from load_and_taylor import load_data_from_mat, fit_polynomial, create_taylor_approximation, plot_approximation, k_lambda
from scipy.optimize import fsolve
import pandas as pd


"""  CIERRA LOS OJOS MAX NO VEAS ESTO   """

neff_df = pd.read_csv('/home/jay/Downloads/NeffP12.csv')
filtered_rows = neff_df[~neff_df.isin([333, 555]).any(axis=1)]
neffLimpio = filtered_rows.copy()
Neff_limpio = neffLimpio.values
lambda_values = np.array(neffLimpio.columns, dtype=float)
Neff_limpio, lambda_values
neff_values = Neff_limpio[0]

def DK(LAMP, LAMS, coeff):
        LAMI = 1. / (2. / LAMP - 1. / LAMS)
        
        kp = k_lambda(LAMP, coeff)
        ks = k_lambda(LAMS, coeff)
        ki = k_lambda(LAMI, coeff)
        
        return kp + kp - ks - ki - 1e-6


def lsli(lamp, NL, as_, coeff):
    lams = fsolve(DK, as_, args=(lamp, coeff), xtol=1e-20, maxfev=10000)[0]
    lami = (2 * 3.14159 * 3e14) / ((2 * 3.14159 * 3e14) / lamp + (2 * 3.14159 * 3e14) / lamp - (2 * 3.14159 * 3e14) / lams)

    return lams, lami

path = '/home/jay/repos/F3001C_Reto/CodigoFinal/Modos/Modes/SweepOverlapTE/Waveguide1000_550_1580.mat'

""" TAYLOR APPROXIMATION """

neff_values, lambda_values = load_data_from_mat(path)

# Fit polynomial
coefficients = fit_polynomial(lambda_values, neff_values)

# Create Taylor approximation
polynomial = np.poly1d(coefficients)
center_lambda = (lambda_values.min() + lambda_values.max()) / 2
taylor_polynomial = create_taylor_approximation(polynomial, center_lambda, order = 30)

# Evaluate polynomials
lambda_range = np.linspace(lambda_values.min(), lambda_values.max(), 500)
original_values = np.polyval(coefficients, lambda_range)
taylor_values = taylor_polynomial(lambda_range - center_lambda)

# Plot Neff approximation
#plot_approximation(lambda_values, neff_values, lambda_range, original_values, taylor_values, 'Neff', 'Taylor Polynomial Approximation of Neff')

# Calculate K(lambda)
K_values = (2 * np.pi * neff_values) / lambda_values
K_original = (2 * np.pi * original_values) / lambda_range
K_taylor = (2 * np.pi * taylor_values) / lambda_range

# Plot K(lambda) approximation
#plot_approximation(lambda_values, K_values, lambda_range, K_original, K_taylor, 'K(Lambda)', 'Taylor Polynomial Approximation of K(Lambda)')


""" PHASE MATHING COMPUTE """

LAMBDA_MIN = 0.3
LAMBDA_MAX = 2.3
POINTS = 3000


# Lambda arrangements
lamp = np.linspace(LAMBDA_MIN, LAMBDA_MAX, POINTS)
lams = np.linspace(LAMBDA_MIN, LAMBDA_MAX, POINTS)

# Lambda Meshgrids
LAMP, LAMS = np.meshgrid(lamp, lams)

# DK_TE = DK(LAMP, LAMS, coefficients)

# # Check the shape and a sample value to ensure calculations are correct
# DK_TE.shape, DK_TE[1000, 1000]

# # Plotting the contour for DK_TE without the colorbar
# plt.figure(figsize=(10, 8))
# plt.contour(lamp, lams, DK_TE, [0], colors='b', linewidths=2)
# plt.title('Contour plot of Phase Mismatch (DK_TE)')
# plt.xlabel('lamp')
# plt.ylabel('lams')
# plt.grid(True)
# plt.tight_layout()
# plt.show()


L = 0.1e6
sigma = 0.1e12
lamp0 = 0.762
NL = 0

Ns = 1000

omgp0 = (2*np.pi*3e14)/lamp0
omgp = np.linspace(omgp0-3*sigma, omgp0+3*sigma,Ns)

a_s = 0.6

[lams0, lami0] = lsli(lamp0, NL, a_s, coefficients)

omgs0 = (2*np.pi*3e14)/lams0
omgi0 = (2*np.pi*3e14)/lami0

dw = 24e12


oms = np.linspace(omgs0 - dw, omgs0 + dw, Ns)
omi = np.linspace(omgi0 - dw, omgi0 + dw, Ns)

[OMS, OMI] = np.meshgrid(oms,omi)

KS = k_lambda((2*np.pi*3e14)/OMS, coefficients)
KI = k_lambda((2*np.pi*3e14)/OMI, coefficients)
KP1 = k_lambda((2*np.pi*3e14)/omgp, coefficients)

dwp = omgp[1]- omgp[0]
JSA = 0

# Integrate JSA
# JSJSJS FOR LOOP TAN RANDOM QUE NI YO LO ENTIENDO

for i in range(0, len(omgp)):
    KP2 = k_lambda(2*np.pi*3e14/(OMS + OMI - omgp[i]), coefficients)
    D_K = KP1 + KP2 - KS - KI - NL 

    # INTEGRAL GO BRRR

    JSA += dwp * np.exp(
                        -((omgp[i]-omgp0)**2)/(sigma**2) *
                        np.exp(-(OMS+OMI-omgp[i]-omgp0)**2 / (sigma**2))
            
                        ) * np.sinc( 
                                    (L/2) * D_K)* np.exp(1j * L * (KP1[i] + KP2 - KS - KI -NL )/2 )
        

plt.figure(figsize=(10, 8))
plt.pcolor(oms, omi, np.abs(JSA)**2)
plt.title('JSA')
plt.xlabel('oms')
plt.ylabel('omi')
plt.grid(True)
plt.tight_layout()
plt.show()


