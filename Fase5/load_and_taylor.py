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

def k_prime_lambda(lambda_val, coefficients):
    """Compute the derivative of the propagation constant for given lambda values using polynomial coefficients."""
    # Compute neff for given lambda_val
    neff = np.polyval(coefficients, lambda_val)
    
    # Derive the coefficients of the polynomial to get the coefficients of neff' (derivative of neff)
    derivative_coefficients = np.polyder(coefficients)
    
    # Evaluate the derivative polynomial to get the value of neff' at lambda_val
    neff_prime = np.polyval(derivative_coefficients, lambda_val)
    
    # Compute k' using the formula above
    k_prime = (2 * np.pi * (neff_prime * lambda_val - neff)) / lambda_val**2
    
    return k_prime

