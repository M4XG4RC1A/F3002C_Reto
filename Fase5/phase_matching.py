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

def pade_coefficients(taylor_coeffs, m, n):
    """
    Compute the Pade approximation coefficients for a given Taylor series.
    
    Parameters:
    - taylor_coeffs: Coefficients of the Taylor series, starting from the 0th degree term.
    - m: Degree of the numerator polynomial P(x).
    - n: Degree of the denominator polynomial Q(x).
    
    Returns:
    - p_coeffs: Coefficients of the numerator polynomial P(x).
    - q_coeffs: Coefficients of the denominator polynomial Q(x).
    """
    
    # Total number of coefficients in the Pade approximation
    N = m + n + 1
    
    # Ensure we have enough Taylor coefficients
    if len(taylor_coeffs) < N:
        raise ValueError("Not enough Taylor coefficients provided.")
    
    # Construct the matrix A and vector B for Ax = B
    L = np.zeros((N, N))
    R = np.zeros(N)

    # First vector L is filled with 0
    for i in range(N):
        L[i] = 0
    # Matrix R is divided in two sections, first section
    # C_n+1 to C_2n
    # B_1 to B_n
    for i in range(N):
        for j in range(N):
            R[i, j] = taylor_coeffs[n+i-j]

    # Solve for B_n
    B = np.linalg.solve(L, R)
    print(B)



    # 
    # 0 to A
    # 0 to C_n  
    # B_1 to B_n




    # and C_n to C_2n


    
    


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

def pade_lambda(lambda_val, p_coeffs, q_coeffs):
    """Compute the Pade approximation for given lambda values using polynomial coefficients."""
    # Evaluate the numerator polynomial P(x) at lambda_val
    p_value = np.polyval(p_coeffs, lambda_val)
    
    # Evaluate the denominator polynomial Q(x) at lambda_val
    q_value = np.polyval(q_coeffs, lambda_val)
    
    # Compute the value of the Pade approximation at lambda_val
    pade_value = p_value / q_value
    
    return pade_value




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
pathTE = '/home/jay/repos/F3001C_Reto/CodigoFinal/Modos/Modes/SweepOverlapTE/Waveguide1000_550_1580.mat'


#pathTM = '/home/jay/repos/F3001C_Reto/CodigoFinal/Modos/Modes/TM/Waveguide1000_325_1555.mat'



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

