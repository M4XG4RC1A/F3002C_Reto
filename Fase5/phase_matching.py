import matplotlib.pyplot as plt
import numpy as np
import h5py
from load_and_taylor import load_data_from_mat, fit_polynomial, create_taylor_approximation, plot_approximation, k_lambda


def compute_k_l(path):

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


    """ PADE APPROXIMATION LETS GOOOO """ 

    #pade_coefficients(coefficients, int(30/2), int(30/2))
    
    # Evaluate polynomial
    #lambda_range = np.linspace(lambda_values.min(), lambda_values.max(), 500)
    #original_values = np.polyval(coefficients, lambda_range)
    #pade_values = pade_lambda(lambda_range, p_coeff, q_coeff)

    # Plot Neff approximation
    #plot_approximation(lambda_values, neff_values, lambda_range, original_values, pade_values, 'Neff', 'Pade Polynomial Approximation of Neff')





    """ PHASE MATHING COMPUTE """

    LAMBDA_MIN = 0.3
    LAMBDA_MAX = 2.3
    POINTS = 3000


    # Lambda arrangements
    lamp = np.linspace(LAMBDA_MIN, LAMBDA_MAX, POINTS)
    lams = np.linspace(LAMBDA_MIN, LAMBDA_MAX, POINTS)

    # Lambda Meshgrids
    LAMP, LAMS = np.meshgrid(lamp, lams)
    LAMI = 1. / (2. / LAMP - 1. / LAMS)

    # Propagation constants for each field
    kp = k_lambda(LAMP, coefficients)
    ks = k_lambda(LAMS, coefficients)
    ki = k_lambda(LAMI, coefficients)

    # # Pade approximation for each field
    # kp = pade_lambda(LAMP, p_coeff, q_coeff)
    # ks = pade_lambda(LAMS, p_coeff, q_coeff)
    # ki = pade_lambda(LAMI, p_coeff, q_coeff)


    return kp, ks, ki, lamp, lams


# bueno
pathTE = '/home/jay/repos/F3001C_Reto/CodigoFinal/Modos/Modes/SweepOverlapTE/Waveguide1000_325_1580.mat'

# buenisimo
pathTE = '/home/jay/repos/F3002C_Reto/Fase4/Sweep/Matlab/Waveguide727778_1000000_1580_Mode3.mat'


#pathTM = '/home/jay/repos/F3001C_Reto/CodigoFinal/Modos/Modes/TM/Waveguide1000_325_1555.mat'

pathTE = 'Fase4/Sweep/Matlab/Waveguide1000000_750000_1580_Mode3.mat'

# Phase mismatch
kp_TE, ks_TE, ki_TE, lamp, lams = compute_k_l(pathTE)
#kp_TM, ks_TM, ki_TM = compute_k_l(pathTM)

DK_TE = kp_TE + kp_TE - ks_TE - ki_TE - 1e-6

# Check the shape and a sample value to ensure calculations are correct
DK_TE.shape, DK_TE[1000, 1000]

# Plotting the contour for DK_TE without the colorbar
plt.figure(figsize=(10, 8))
plt.contour(lamp, lams, DK_TE, [0], colors='b', linewidths=2)
plt.title('Contour plot of Phase Mismatch (DK_TE)')
plt.xlabel('lamp')
plt.ylabel('lams')
plt.grid(True)
plt.tight_layout()
plt.show()
