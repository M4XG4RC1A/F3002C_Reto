import numpy as np
import matplotlib.pyplot as plt
import math

def fit_polynomial(lambda_values, neff_values, degree=30):
    """Fit a polynomial of given degree and return coefficients."""
    return np.polyfit(lambda_values, neff_values, degree)


def create_taylor_approximation(polynomial, center, order=30):
    """Create a Taylor polynomial approximation centered around a given value."""
    taylor_coefficients = [polynomial(center)]
    for i in range(1, order+1):
        polynomial = np.polyder(polynomial)
        taylor_coefficient = polynomial(center) / math.factorial(i)
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
    N = m + n
    
    # Ensure we have enough Taylor coefficients
    if len(taylor_coeffs) < N:
        raise ValueError("Not enough Taylor coefficients provided.")
    
    # Construct the matrix A and vector B for Ax = B
    R = np.zeros((N, N))
    L = np.zeros(N)
    

    """     
    THE MATRIX TO SOLVE IS DIVIDED INTO 4 SECTIONS, m-1 B_n   ->   n+1 A_n = 0
                                    2nd equation    C_n B_n   ->   -A_n
    """

    # First stage, L is filled with -C_n+1 to -C_2n
    for i in range(m):
        L[i] = -taylor_coeffs[m+i+1]

    # Next stage, other half of L is filled with C_n
    for i in range(n):
        L[i+m] = -taylor_coeffs[i+1]

    # First stage fill top left matrix
    # C_n+1 to C_2n
    # B_1 to B_n
    for i in range(m):
        for j in range(m):
            R[i,j] = taylor_coeffs[n+i-j]
    

    # TOP RIGHT IS FILLED WITH 0 BECAUSE THERE ARE NO A_N

    # BOTTOM LEFT
    for i in range(n):
        for j in range(m):
            if i-j >= 0:
                R[i+n,j] = taylor_coeffs[i-j]

    # LAST STEP WILL WITH -A_N
    for i in range(n):
        R[m+i, m+i] = -1

    # # Solve for B_n
    Solutions = np.linalg.solve(R, L)
    #print(taylor_coeffs[0])

    A = np.concatenate(([taylor_coeffs[0]], Solutions[n:]))
    B = Solutions[:n]

    return A, B
    
def eval_pade(x, p_coeffs, q_coeffs):
    denominator = 1
    for i in range(len(q_coeffs)):
        denominator += q_coeffs[i] * x**i
    numerator = 0
    for i in range(len(p_coeffs)):
        numerator += p_coeffs[i] * x**i
    return numerator / denominator


f = lambda x : x**2 + 2*x + 1 + np.random.normal(0, 0.4, len(x))

x = np.linspace(0, 5, 10)
y = f(x)


order = 30

coefficients = fit_polynomial(x, y, degree=order)
polynomial = np.poly1d(coefficients)

taylor = create_taylor_approximation(polynomial, 2, order=order)
taylor_values = taylor(x-2)

p_coefficients = pade_coefficients(coefficients, 1, 1)


# print(p_coefficients)
pade = eval_pade(x, p_coefficients[0], p_coefficients[1])




plt.scatter(x,y)
plt.plot(x, pade)
plt.show()

