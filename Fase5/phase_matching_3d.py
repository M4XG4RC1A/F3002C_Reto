import matplotlib.pyplot as plt
import numpy as np
import h5py
from load_and_taylor import load_data_from_mat, fit_polynomial, create_taylor_approximation, plot_approximation, k_lambda


def compute_k_l(path, two_D=False):
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

    #  # Calculate K(lambda)
    # K_values = (2 * np.pi * neff_values) / lambda_values
    # K_original = (2 * np.pi * original_values) / lambda_range
    # K_taylor = (2 * np.pi * taylor_values) / lambda_range

    # # Plot K(lambda) approximation
    # plot_approximation(lambda_values, K_values, lambda_range, K_original, K_taylor, 'K(Lambda)', 'Taylor Polynomial Approximation of K(Lambda)')


    """

    Begin 3d arrangement of DK 
    
    """
    LAMBDA_MIN = 0
    LAMBDA_MAX = 2
    POINTS = 300

    # Lambda arrangements
    lamp = np.linspace(LAMBDA_MIN, LAMBDA_MAX, POINTS)
    lams = np.linspace(LAMBDA_MIN, LAMBDA_MAX, POINTS)
    lamr = np.linspace(LAMBDA_MIN, LAMBDA_MAX, POINTS)

    # Lambda Meshgrids
    LAMP, LAMS, LAMR = np.meshgrid(lamp, lams, lamr)

    LAMI = 1. / (1. / LAMP - 1. / LAMS - 1. / LAMR)

    # Propagation constants for each field
    kp = k_lambda(LAMP, coefficients)
    ks = k_lambda(LAMS, coefficients)
    kr = k_lambda(LAMR, coefficients)
    ki = k_lambda(LAMI, coefficients)


    

    return kp, ks, kr, ki, lamp, lams, lamr


def hex_to_rgb(value):
    """Convert hex string to a tuple of RGB."""
    value = value.lstrip('#')
    length = len(value)
    return tuple(int(value[i:i+length//3], 16) for i in range(0, length, length//3))

def rgb_to_hex(rgb):
    """Convert a tuple of RGB values to a hex string."""
    return '#{:02x}{:02x}{:02x}'.format(*rgb)

def interpolate_color(color1, color2, factor):
    """Interpolate between two RGB colors."""
    r1, g1, b1 = color1
    r2, g2, b2 = color2
    r = r1 + (r2 - r1) * factor
    g = g1 + (g2 - g1) * factor
    b = b1 + (b2 - b1) * factor
    return int(r), int(g), int(b)



# # bueno
# pathTE = '/home/jay/repos/F3001C_Reto/CodigoFinal/Modos/Modes/SweepOverlapTE/Waveguide1000_325_1580.mat'

# # buenisimo
# pathTE = '/home/jay/repos/F3002C_Reto/Fase4/Sweep/Matlab/Waveguide727778_1000000_1580_Mode3.mat'


pathTE = '/home/jay/repos/F3002C_Reto/Fase4/Sweep/Matlab/Waveguide727778_1000000_1580_Mode3.mat'

pathTE = '/home/jay/repos/F3002C_Reto/Fase4/Sweep/Matlab/Waveguide727778_1000000_1580_Mode3.mat'




#pathTM = '/home/jay/repos/F3001C_Reto/CodigoFinal/Modos/Modes/TM/Waveguide1000_325_1555.mat'

# Phase mismatch
kp, ks, ki, kr, lamp, lams, lamr = compute_k_l(pathTE)
#kp_TM, ks_TM, ki_TM = compute_k_l(pathTM)

DK_2d = kp + kp - ks - ki - 1e-6

DK_3d = kp - ks - kr - ki - 1e-6

# Plotting the contour for DK_TE without the colorbar
plt.figure(figsize=(10, 8))

start_color = hex_to_rgb("#FF0000")
end_color = hex_to_rgb("#0000FF")

n = 100
for i in range(n):
    slice_index = int(len(DK_3d)*i/n)
    
    # Get interpolated color
    factor = i/n
    rgb_color = interpolate_color(start_color, end_color, factor)
    hex_color = rgb_to_hex(rgb_color)

    # Here's your plotting line with the new color:
    plt.contour(lamp, lams, DK_3d[:, slice_index, :], 0, colors=hex_color, linewidths=0.5)


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

