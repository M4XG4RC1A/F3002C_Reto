import matplotlib.pyplot as plt
import numpy as np
import h5py

def load_data_from_mat(file_path):
    """Load data from .mat file and return neff and lambda values."""
    with h5py.File(file_path, 'r') as file:
        neff_values = np.array(file['neff'][0])
        lambda_values = np.array(file['lambda'][0])
    return neff_values, lambda_values

def fit_polynomial(lambda_values, neff_values, degree=30):
    """Fit a polynomial of given degree and return coefficients."""
    return np.polyfit(lambda_values, neff_values, degree)

def create_taylor_approximation(polynomial, center, order=30):
    """Create a Taylor polynomial approximation centered around a given value."""
    taylor_coefficients = [polynomial(center)]
    for i in range(1, order+1):
        polynomial = np.polyder(polynomial)
        taylor_coefficient = polynomial(center) / np.math.factorial(i)
        taylor_coefficients.append(taylor_coefficient)
    return np.poly1d(taylor_coefficients[::-1])


def plot_approximation(lambda_values, neff_values, lambda_range, original_values, taylor_values, ylabel, title):
    """Plot the original data, fitted polynomial, and Taylor approximation."""
    plt.figure(figsize=(12, 6))
    plt.plot(lambda_values, neff_values, 'o', label='Original Data')
    plt.plot(lambda_range, original_values, '--', label='Fitted Polynomial')
    plt.plot(lambda_range, taylor_values, '-', label='Taylor Approximation')
    plt.xlabel('Lambda (um)')
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def k_lambda(lambda_val, coefficients):
    """Compute the propagation constant for given lambda values using polynomial coefficients."""
    neff = np.polyval(coefficients, lambda_val)
    return (2 * np.pi * neff) / lambda_val


def compute_k_l(path):
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

    # Calculate K(lambda)
    K_values = (2 * np.pi * neff_values) / lambda_values
    K_original = (2 * np.pi * original_values) / lambda_range
    K_taylor = (2 * np.pi * taylor_values) / lambda_range

    # Plot K(lambda) approximation
    plot_approximation(lambda_values, K_values, lambda_range, K_original, K_taylor, 'K(Lambda)', 'Taylor Polynomial Approximation of K(Lambda)')


    """

    Begin 3d arrangement of DK 
    
    """


    # Lambda arrangements
    lamp = np.linspace(0.7, 1.7, 30)
    lams = np.linspace(0.7, 1.7, 30)
    lamr = np.linspace(0.7, 1.7, 30)

    # Lambda Meshgrids
    LAMP, LAMS, LAMR = np.meshgrid(lamp, lams, lamr)

    LAMI = 1. / (2. / LAMP - 1. / LAMS - 1. / LAMR)

    # Propagation constants for each field
    kp = k_lambda(LAMP, coefficients)
    ks = k_lambda(LAMS, coefficients)
    kr = k_lambda(LAMR, coefficients)
    ki = k_lambda(LAMI, coefficients)

    return kp, ks, kr, ki, lamp, lams, lamr


# bueno
pathTE = '/home/jay/repos/F3001C_Reto/CodigoFinal/Modos/Modes/SweepOverlapTE/Waveguide1000_325_1580.mat'

# buenisimo
pathTE = '/home/jay/repos/F3001C_Reto/CodigoFinal/Modos/Modes/SweepOverlapTE/Waveguide1000_550_1580.mat'


#pathTM = '/home/jay/repos/F3001C_Reto/CodigoFinal/Modos/Modes/TM/Waveguide1000_325_1555.mat'



# Phase mismatch
kp, ks, ki, kr, lamp, lams, lamr = compute_k_l(pathTE)
#kp_TM, ks_TM, ki_TM = compute_k_l(pathTM)

DK_TE = kp + kp - ks - kr - ki - 1e-6

# # Check the shape and a sample value to ensure calculations are correct
# DK_TE.shape, DK_TE[30, 30, 30]

# Plotting the contour for DK_TE without the colorbar
plt.figure(figsize=(10, 8))

# Plot 3D surface
plt.contour(lamp, lams, DK_TE[:, 15, :], 100, cmap='jet')


plt.title('Contour plot of Phase Mismatch (DK_TE)')
plt.xlabel('lamp')
plt.ylabel('lams')
plt.grid(True)
plt.tight_layout()
plt.show()



# # Phase mismatch
# kp, ks, ki, kr, lamp, lams, lamr = compute_k_l(pathTE)
# #kp_TM, ks_TM, ki_TM = compute_k_l(pathTM)

# DK_TE = kp + kp - ks - kr - ki - 1e-6

# # # Check the shape and a sample value to ensure calculations are correct
# # DK_TE.shape, DK_TE[30, 30, 30]

# # Plotting the contour for DK_TE without the colorbar
# fig = plt.figure(figsize=(10, 8))
# ax = fig.add_subplot(111, projection="3d")
# ax.plot_surface(lamr, lams, DK_TE, cmap="autumn_r", lw=0, rstride=1, cstride=1)

# ax.title('Contour plot of Phase Mismatch (DK_TE)')
# ax.xlabel('lamp')
# ax.ylabel('lams')
# ax.grid(True)
# ax.tight_layout()
# plt.show()

