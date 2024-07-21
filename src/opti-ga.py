import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from numpy.polynomial.legendre import Legendre
from IPython.display import clear_output
import time

# NUCLEAR PARAMETERS (same as in your original code)
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
M = 6
dx = L / (N-1)
x = np.linspace(0, L, N)

# GA-RELATED PARAMETERS
lower_bound = -0.1
upper_bound = 0.1
epsilon = 0.1
pop_size=10
num_generations=50
mutation_rate = 0.4
sigma = 0.01
alpha = 0.3

## COST-FUNCTION LAMBDAS
lambda_ppf = 1e3
lambda_flux = 1e1
lambda_avg = 1e1

def map_to_legendre_interval(x, L):
    """Map x from [0, L] to [-1, 1]."""
    return 2 * x / L - 1

def legendre_polynomials_on_interval(x, n, L):
    """
    Evaluate the first n Legendre polynomials on the interval [0, L].
    
    Parameters:
    x (array_like): Points at which to evaluate the polynomials.
    n (int): Number of Legendre polynomials to evaluate.
    L (float): The upper bound of the interval [0, L].
    
    Returns:
    list of arrays: Each array contains the values of the nth Legendre polynomial at x.
    """
    xi = map_to_legendre_interval(x, L)
    polys = [Legendre.basis(i)(xi) for i in range(n)]
    return polys

def linear_combination_legendre(x, coefficients, L):
    """
    Compute a linear combination of Legendre polynomials on [0, L].
    
    Parameters:
    x (array_like): Points at which to evaluate the polynomials.
    coefficients (array_like): Coefficients for the linear combination.
    L (float): The upper bound of the interval [0, L].
    
    Returns:
    array: Values of the linear combination at x.
    """
    n = len(coefficients)
    polynomials = legendre_polynomials_on_interval(x, n, L)
    linear_combination = sum(c * p for c, p in zip(coefficients, polynomials))
    return linear_combination

def construct_matrix_A(perc):
    """
    Constructs the matrix A for the neutron diffusion equation.

    Args:
        perc (float): Percentage of control rod insertion.

    Returns:
        numpy.ndarray: Constructed matrix A.
    """
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
    """
    Performs the inverse power iteration to find the dominant eigenvalue and eigenvector.

    Args:
        A (numpy.ndarray): Matrix A.
        B (numpy.ndarray): Matrix B.
        tol (float): Convergence tolerance. Defaults to 1e-4.

    Returns:
        float: Dominant eigenvalue.
        numpy.ndarray: Corresponding eigenvector.
    """
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
    """
    Solves for the zeroth-order neutron flux.

    Args:
        perc (float): Percentage of control rod insertion.

    Returns:
        numpy.ndarray: Zeroth-order neutron flux.
    """
    A = construct_matrix_A(perc)
    b = nu * Sigma_f * np.ones((N-1, N-1))
    k, phi = invPow(A, b)
    phi[0] = phi[-1] = 0
    return phi

def solve_phi1(N1, perc):
    """
    Solves for the first-order perturbation solution for neutron flux.

    Args:
        N1 (numpy.ndarray): Initial perturbation function.
        perc (float): Percentage of control rod insertion.

    Returns:
        numpy.ndarray: First-order perturbation solution for neutron flux.
    """

    ## implement n:th order perturbation
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
    """
    Objective function to minimize.

    Args:
        c_nu (numpy.ndarray): Coefficients for Legendre polynomials.

    Returns:
        float: Calculated ppf loss, flux loss, avg loss and total loss.
    """
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

    for perc in np.arange(0.0, 1, 0.1):
        phi0 = solve_phi0(perc)
        phi1 = solve_phi1(N1, perc)
        phi = phi0 + epsilon * phi1
        ppf0 = np.max(phi0) / np.mean(phi0)
        ppf = np.max(phi) / np.mean(phi)
        ppf_0s.append(ppf0)
        ppfs.append(ppf)
        loss_ppf = (ppf - ppf0)
        avg = abs(np.sum((N1) * dx))
        flux = abs(np.sum((phi-phi0) * dx))
        ppf_loss += loss_ppf * lambda_ppf
        flux_loss += flux * lambda_flux
        avg_loss += avg * lambda_avg

        loss += loss_ppf * lambda_ppf + lambda_flux * flux + lambda_avg * avg

    return ppf_loss, flux_loss, avg_loss, loss 

def initialize_population(pop_size, n):
    """ Initialize population of size pop_size with n-dimensional individuals."""
    return np.random.uniform(lower_bound, upper_bound, size=(pop_size, n))

def select_parents(population, fitness):
    """ Select parents based on fitness (lower fitness is better)."""
    sorted_indices = np.argsort(fitness)
    return population[sorted_indices[:2]]

def crossover(parent1, parent2, alpha=alpha):
    """ Perform crossover between two parents to produce offspring."""
    child1 = alpha * parent1 + (1 - alpha) * parent2
    child2 = alpha * parent2 + (1 - alpha) * parent1
    return child1, child2

def mutate(individual, mutation_rate=mutation_rate, sigma=sigma):
    """
    Apply mutation to an individual in the population.

    Args:
        individual (numpy.ndarray): The individual to be mutated.
        mutation_rate (float): Probability of mutating each gene in the individual.
        sigma (float): Initial standard deviation for the mutation noise.

    Returns:
        numpy.ndarray: The mutated individual.
    """
    for i in range(len(individual)):
        if np.random.rand() < mutation_rate:
            individual[i] += np.random.normal(0, sigma)
    return individual

def ga_optimize(pop_size=pop_size, num_generations=num_generations):

    """
    Perform genetic algorithm optimization to find the best solution for the control rod design.

    Args:
        pop_size (int): The size of the population.
        num_generations (int): The number of generations to run the optimization.

    Returns:
        tuple: A tuple containing:
            - best_solution (numpy.ndarray): The best solution found.
            - best_fitness (float): The fitness of the best solution.
            - fitness_over_time (list of floats): The fitness values over generations.
    """
     
    population = initialize_population(pop_size, M)
    best_solution = None
    best_fitness = float('inf')

    fitness_over_time = []
    ppf_over_time = []
    avg_over_time = []
    flux_over_time = []

    plt.ion()
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_xlabel('Generation', size = '16')
    ax.set_ylabel(r"$\log \mathcal{C}(\mathbf{c_\eta})$", size = '16')
    ax.grid(True)

    ppf_line, = ax.plot([], [], 'r-', label='PPF Fitness')
    avg_line, = ax.plot([], [], 'y-', label='Average Fitness')
    flux_line, = ax.plot([], [], 'g-', label='Flux Fitness')
    total_line, = ax.plot([], [], 'm-', label='Total Fitness')
    ax.set_xlim(0, num_generations)
    ax.legend()

    progress_bar = tqdm(range(num_generations), desc="Optimizing", leave=True)

    for generation in progress_bar:
        ppf_fitness, flux_fitness, avg_fitness, total_fitness = zip(*[objective(individual) for individual in population])
        fitness = np.array([objective(individual)[-1] for individual in population])  # Total fitness
        ppf_fitness = np.array([objective(individual)[0] for individual in population])
        flux_fitness = np.array([objective(individual)[1] for individual in population])
        avg_fitness = np.array([objective(individual)[2] for individual in population])

        ppf_over_time.append(np.min(ppf_fitness))
        flux_over_time.append(np.min(flux_fitness))
        avg_over_time.append(np.min(avg_fitness))
        fitness_over_time.append(np.min(fitness))

        parent1, parent2 = select_parents(population, fitness)

        new_population = []
        for _ in range(pop_size // 2):
            child1, child2 = crossover(parent1, parent2)
            child1 = mutate(child1)
            child2 = mutate(child2)
            new_population.extend([child1, child2])

        population = np.array(new_population)

        for _ in range(pop_size // 2):
            child1, child2 = crossover(parent1, parent2)
            child1 = mutate(child1)
            child2 = mutate(child2)
            new_population.extend([child1, child2])

        population = np.array(new_population)

        if np.min(fitness) < best_fitness:
            best_fitness = np.min(fitness)
            best_solution = population[np.argmin(fitness)]

        progress_bar.set_postfix(best_fitness=best_fitness)
        progress_bar.set_postfix(best_fitness=best_fitness)

        # Update live plot
        ppf_line.set_xdata(range(len(fitness_over_time)))
        ppf_line.set_ydata(ppf_over_time)
        flux_line.set_xdata(range(len(fitness_over_time)))
        flux_line.set_ydata(flux_over_time)
        avg_line.set_xdata(range(len(fitness_over_time)))
        avg_line.set_ydata(avg_over_time)
        total_line.set_xdata(range(len(fitness_over_time)))
        total_line.set_ydata(fitness_over_time)

        ax.set_ylim(np.min([min(ppf_over_time), min(fitness_over_time),  min(flux_over_time), min(avg_over_time)]),
                     np.max([max(ppf_over_time), max(fitness_over_time),  max(flux_over_time), max(avg_over_time)])) 
        fig.canvas.draw()
        fig.canvas.flush_events()
        ax.set_xlim(0, generation)
        clear_output(wait=True)
        display(fig)
        time.sleep(0.01)

    plt.ioff()
    return best_solution, best_fitness, fitness_over_time
     
# Running the Genetic Algorithm optimization
best_solution, best_fitness, fitness_over_time = ga_optimize()
