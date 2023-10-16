import matplotlib.pyplot as plt
import numpy as np
import h5py
from load_and_taylor import load_data_from_mat, fit_polynomial, create_taylor_approximation, plot_approximation, k_lambda
from scipy.optimize import fsolve
import pandas as pd
import os

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

def plot_contour(lamp, lams, DK_TE, path):

    plt.figure(figsize=(10, 8))
    plt.contour(lamp, lams, DK_TE, [0], colors='b', linewidths=2)
    plt.title('Contour plot of Phase Mismatch (DK_TE)')
    plt.xlabel('lamp')
    plt.ylabel('lams')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.dirname(path) + '/DK' + os.path.basename(path) + '.png')

def plot_correlation(oms, omi, F_cp, path):

    plt.figure(figsize=(10, 8))
    plt.pcolor(oms, omi, np.abs(F_cp)**2)
    # Add colorbar
    cbar = plt.colorbar()
    cbar.set_label('Intensity')
    plt.title('JSA')
    plt.xlabel('oms')
    plt.ylabel('omi')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.dirname(path) + '/F_cp' + os.path.basename(path) + '.png')


def correlation(path):

    
    """ TAYLOR APPROXIMATION """

    neff_values, lambda_values = load_data_from_mat(path)

    """  CIERRA LOS OJOS MAX NO VEAS ESTO   """

    neff_df = pd.read_csv('/home/jay/Downloads/NeffP12.csv')
    filtered_rows = neff_df[~neff_df.isin([333, 555]).any(axis=1)]
    neffLimpio = filtered_rows.copy()
    Neff_limpio = neffLimpio.values
    lambda_values = np.array(neffLimpio.columns, dtype=float)
    Neff_limpio, lambda_values
    neff_values = Neff_limpio[0]
    
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

    DK_TE = DK(LAMP, LAMS, coefficients)

    # # # Check the shape and a sample value to ensure calculations are correct
    DK_TE.shape, DK_TE[1000, 1000]

    plot_contour(lamp, lams, DK_TE, path)

    L = 0.1e6 # in mm
    SIGMA = 0.1e12
    L0 = 0.798
    NL = 0

    N = 500

    omp0 = (2*np.pi*3e14)/L0
    omp = np.linspace(omp0-3*SIGMA, omp0+3*SIGMA,N)

    a_s = 0.6

    [lams0, lami0] = lsli(L0, NL, a_s, coefficients)

    oms0 = (2*np.pi*3e14)/lams0
    omi0 = (2*np.pi*3e14)/lami0

    dw = 24e12


    oms = np.linspace(oms0 - dw, oms0 + dw, N)
    omi = np.linspace(omi0 - dw, omi0 + dw, N)

    [OMS, OMI] = np.meshgrid(oms,omi)

    KS = k_lambda((2*np.pi*3e14)/OMS, coefficients)
    KI = k_lambda((2*np.pi*3e14)/OMI, coefficients)
    KP1 = k_lambda((2*np.pi*3e14)/omp, coefficients)

    d_wp = omp[1]- omp[0]
    F_cp = 0

    def alpha1 (omg, sig):
        return -((omg)**2)/(sig**2)

    # integral go brrrrrr
    for j in range(len(omp)):
        KP2 = k_lambda(2*np.pi*3e14/(OMS + OMI - omp[j]), coefficients)
        delta_K = KP1[j] + KP2 - KS - KI - NL

        F_cp += d_wp *(
                    np.exp(alpha1(omp[j] - omp0,SIGMA)) *
                    np.exp(alpha1(OMS + OMI - omp[j] - omp0,SIGMA)) *
                    np.sinc(L * (delta_K) / 2)  *
                    np.exp(1j * L * (delta_K) / 2)
                    )
        
    plot_correlation(oms,omi,F_cp,path)


mode = '/home/jay/repos/F3002C_Reto/Fase4/Sweep/Matlab/'
mode = '/home/jay/repos/F3002C_Reto/Fase3/Sweeps/Matlab/'

# Loop over all files inside mode directory
for file in os.listdir(mode):
    correlation(mode + file)

