import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from numpy.polynomial.legendre import Legendre

# NUCLEAR PARAMETERS
D = 3.85
Sigma_a = 0.0532
nu = 2.5
k = 1.0
Sigma_f = 0.157
sigma_c = 200

# SPATIAL PARAMETERS
N0 = 1
L = 1.0
N = 50
M = 5
dx = L / (N-1)
x = np.linspace(0, L, N)

# GA-RELATED PARAMETERS
lower_bound = -1
upper_bound = 1
epsilon = 0.1

## COST-FUNCTION LAMBDAS
lambda_ppf = 1e1
lambda_flux = 1e1
lambda_avg = 1

def map_to_legendre_interval(x, L):
    return 2 * x / L - 1

def legendre_polynomials_on_interval(x, n, L):
    xi = map_to_legendre_interval(x, L)
    polys = [Legendre.basis(i)(xi) for i in range(n)]
    return polys

def linear_combination_legendre(x, coefficients, L):
    n = len(coefficients)
    polynomials = legendre_polynomials_on_interval(x, n, L)
    linear_combination = sum(c * p for c, p in zip(coefficients, polynomials))
    return linear_combination

def construct_matrix_A(perc):
    A = np.zeros((N-1, N-1))
    common_factor = 2 * D / dx**2
    for i in range(1, N):
        control_rod_effect = N0 * sigma_c if i / N < perc else 0
        if i == 1:
            A[0][0] = common_factor + Sigma_a + control_rod_effect
            A[0][1] = -common_factor / 2
        elif i == N-1:
            A[N-2][N-3] = -common_factor / 2
            A[N-2][N-2] = common_factor + Sigma_a + control_rod_effect
        else:
            A[i-1][i-2] = -common_factor / 2
            A[i-1][i-1] = common_factor + Sigma_a + control_rod_effect
            A[i-1][i] = -common_factor / 2
    return A

def invPow(A, B, tol=1e-4):
    phi = np.random.random((A.shape[0]))  # initial guess
    C = np.linalg.inv(A) @ B
    converged = False
    kold = 0.0
    while not converged:
        phi = np.dot(C, phi)
        k = np.max(np.abs(phi))
        phi = phi / np.max(phi)
        if abs(kold - k) < tol:
            converged = True
        kold = k
    return k, phi

def solve_phi0(perc):
    A = construct_matrix_A(perc)
    b = nu * Sigma_f * np.ones((N-1, N-1))
    k, phi = invPow(A, b)
    phi[0] = phi[-1] = 0
    return phi

def solve_phi1(N1, perc):
    threshold = int(perc * N)
    S = nu * Sigma_f * np.zeros(N-1)
    phi0 = solve_phi0(perc)
    if threshold > 0:
        S[:threshold] += sigma_c * N1[:threshold] * phi0[:threshold]
    A = construct_matrix_A(perc)
    B = nu * Sigma_f / k * np.ones((N-1, N-1))
    phi1 = -np.linalg.inv(A - B) @ S
    phi1[0] = phi1[-1] = 0
    return phi1

def objective(c_nu):
    loss = 0

    c_nu = np.reshape(c_nu, (M,))
    c_nu = c_nu / np.sqrt(np.sum(c_nu**2))

    N1 = linear_combination_legendre(x, c_nu, L)
    N1 = np.reshape(N1, (N,))

    ppf_loss = 0
    avg_loss = 0
    flux_loss = 0

    ppf_0s = []
    ppfs = []

    saf_margin = 1.8
    penalty = 1e6

    for perc in np.arange(0.0, 1, 0.05):
        phi0 = solve_phi0(perc)
        phi1 = solve_phi1(N1, perc)
        phi = phi0 + epsilon * phi1
        ppf0 = np.max(phi0) / np.mean(phi0)
        ppf = np.max(phi) / np.mean(phi)
        ppf_0s.append(ppf0)
        ppfs.append(ppf)
        loss_ppf = (ppf - ppf0)

        if ppf > saf_margin:
            loss += penalty
        avg = abs(np.sum((N1) * dx))
        flux = abs(np.sum((phi-phi0) * dx))
        ppf_loss += loss_ppf * lambda_ppf
        flux_loss += flux * lambda_flux
        avg_loss += avg * lambda_avg

        loss += loss_ppf * lambda_ppf + lambda_flux * flux + lambda_avg * avg

    return loss

def nelder_mead_optimize():
    initial_guess = np.random.uniform(lower_bound, upper_bound, M)
    #normalize initial guess
    initial_guess = initial_guess / np.linalg.norm(initial_guess)
    result = minimize(objective, initial_guess, method='Nelder-Mead', options={'maxiter': 500, 'disp': True})
    return result.x, result.fun

# Running the Nelder-Mead optimization
best_solution= None
best_fitness = 1e12
for i in range(0,40):
    best_solution_i, best_fitness_i = nelder_mead_optimize()
    if best_fitness_i < best_fitness:
        best_solution = best_solution_i
        best_fitness = best_fitness_i

# Print results
print("Best solution:", best_solution)
print("Best fitness:", best_fitness)
