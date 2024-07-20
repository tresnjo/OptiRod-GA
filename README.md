# OptiRod-GA
An attempt at developing an optimized control rod design for a nuclear reactor using a genetic algorithm and perturbation theory of the simple neutron diffusion equation model.
# Table of Contents
1. [Objective](#objective)
2. [Background](#background)
   1. [Diffusion Equation](#diffusion-equation)
   2. [Modeling the Control Rod](#modeling-the-control-rod)
   3. [Perturbation Theory](#perturbation-theory)
3. [Method](#method)
   1. [Steps in Solution Procedure and Cost Function](#steps-in-solution-procedure-and-cost-function)
   2. [Inferring a Control Rod Design](#inferring-a-control-rod-design)
5. [Code Description](#code-description)
6. [Results](#results)
7. [Conclusion](#conclusion)

   
## Objective
The objective of this project is to study whether we can optimize the control rod design by minimizing the axial power peaking factor throughout the insertion process. We will employ perturbation theory to a simple one dimensional neutron diffusion equation, and infer a control rod design from a mimization problem using a genetic algorithm.

## Background

### Diffusion Equation

To achieve the objective, we will as previously stated, employ a simple diffusion equation for a multiplying medium given by the equation

$$-D\frac{d^2\phi}{dx^2}+\Sigma_a\phi(x)=\frac{\nu}{k}\Sigma_f\phi(x).$$

Here $\phi(x)$ represents the neutron flux, $D$ is the diffusion coefficient, $\Sigma_a$ is the absorption cross-section, $\nu$ is the neutrons created per fission, $k$ is the effective multiplication factor, and $\Sigma_f$ is the fission cross-section. The first term on the LHS represents the spreds of neutrons (diffusivity), determined by the diffusion coefficient. The second term represents losses due to absorption in the medium. The first term on the RHS represents gains due to fission proccesses. 

### Modeling the Control Rod

We can model the existance of a control rod by introducing a new absorption cross-section for the control rod $\Sigma_c$. We can also relate the macroscopic cross section with the microscopic cross section through the number density of absorbing particles as

$$\Sigma_c = N \sigma_c.$$

If we let $N$ to be space dependent along the axial direction $x$, we can model the new diffusion equation for a multiplying medium with a fully inserted control rod by

$$-D\frac{d^2\phi}{dx^2}+(\Sigma_a+N(x)\sigma_c)\phi(x)=\frac{\nu}{k}\Sigma_f\phi(x).$$

However, if we only insert the control rod partially (to say $x_c$), we can refine our model using a Heaviside function $\mathcal{H}$ 

$$-D\frac{d^2\phi}{dx^2}+(\Sigma_a+N(x)(1-\mathcal{H}(x-x_c))\sigma_c)\phi(x)=\frac{\nu}{k}\Sigma_f\phi(x).$$

The idea is to find an $N(x)$, and thereby inferring a control rod design, by optimizing the axial power peaking factor, which is defined by the ratio between peak flux and mean flux, as expressed below

$$\text{PPF} = \frac{\text{max}_{x\in[0,L]} \phi(x)}{\frac{1}{L} \int_0^L \phi(x) dx}.$$

If we opt to minimize this for all insertion distances $x_c$, we can sum over of all the PPF's above. This will be more clear when we define our cost function later.

### Perturbation Theory

So how can we find an optimal $N(x)$? The first thing to realize is that we can solve the differential equation above analytically for $N(x) = N_0$, meaning a constant number density, given some boundary conditions ($\phi(0)=\phi(L) = 0$ in this case). A constant number density $N_0$ could be viewed as inserting a control rod with a constant radius. 

The second thing to realize is that we can use perturbation theory to get a solution for a spatially dependent $N(x) $ by letting $N(x) = N_0 + \sum_{j} \epsilon^j N_j(x)$ where $N_j(x)$ is our perturbed number density of order $j$. $\epsilon$ is a small parameter. In order to use perturbation theory, we must also make sure that $|\sum_j \epsilon^j N_j(x) | << N_0 \forall x $. Similiarly we let $\phi(x) = \phi_0 + \sum_j \epsilon^j \phi_j(x)$ where $\phi_j(x)$ is the perturbed neutron flux of order $j$. 

In what follows, we will only care about the first order perturbation $j = 1$ since the contribution perturbations that follow diminish quickly as $\mathcal{O}(\epsilon^j)$. Inserting the above into the diffusion equation yields

$$-D\frac{d^2(\phi_0(x)+\epsilon\phi_1(x))}{dx^2}+\[ \Sigma_a+ \{N_0+\epsilon N_1(x)\} (1-\mathcal{H}(x-x_c))\sigma_c \] (\phi_0(x)+\epsilon\phi_1(x))=\frac{\nu}{k}\Sigma_f(\phi_0(x)+\epsilon\phi_1(x)).$$

By omitting contributions of order $\mathcal{\epsilon^2}$, we get two seperate differential equations by collecting terms of equal perturbation order. For the zeroth and first order perturbation $\epsilon^0$ and $\epsilon^1$ we therefore have

$$
\begin{align*}
\underline{\epsilon^0}: & \qquad  \qquad -D\frac{d^2\phi_0}{dx^2}+[\Sigma_a+N_0\{1-\mathcal{H}(x-x_c)\}\sigma_c]\phi_0 &=\frac{\nu}{k}\Sigma_f\phi_0(x) \newline
\underline{\epsilon^1}: &\qquad  \qquad -D\frac{d^2\phi_1}{dx^2}+[\Sigma_a+N_0\{1-\mathcal{H}(x-x_c)\}\sigma_c]\phi_1 + \textcolor{red}{N_1(x)(1-\mathcal{H}(x-x_c))\sigma_c\phi_0}&=\frac{\nu}{k}\Sigma_f\phi_1(x).
\end{align*}
$$

We notice that the equations are almost identical in $\phi_0$ and $\phi_1$, except a new combinaton term consisting of $N_1$ and the zeroth order solution $\phi_0$ arises (colored in red). For easier reading, we will define the axial cross section

$$ \sigma(x, x_c) = [1-\mathcal{H}(x-x_c)]\sigma_c $$

resulting in the following simplified equations

$$
\begin{align*}
\underline{\epsilon^0}: & \qquad  \qquad -D\frac{d^2\phi_0}{dx^2}+\[ \Sigma_a + N_0\sigma(x, x_c) \]\phi_0 &=\frac{\nu}{k}\Sigma_f\phi_0(x) \newline
\underline{\epsilon^1}: &\qquad  \qquad -D\frac{d^2\phi_1}{dx^2}+\[ \Sigma_a + N_0\sigma(x, x_c) \]\phi_1 + \textcolor{red}{N_1(x)\sigma(x,x_c)\phi_0}&=\frac{\nu}{k}\Sigma_f\phi_1(x).
\end{align*}
$$

We will furthermore expand $N_1$ in an orthonormal series given by the Legendre polynomials $P_n(\frac{2x}{L}-1)$ on $[0,L]$, meaning

$$N_1(x) = \sum_{\eta = 0}^M c_{\eta} P_\eta(\frac{2x}{L}-1)$$

where we will restrict ourselves to some arbitrary large $M$ and make sure $c_\eta$ are properly normalized. Of course if we choose a larger $M$, we can search through a larger subspace of all possible control rod designs, but at the cost of slower computation time. 

## Method

###  Steps in Solution Procedure and Cost Function
The solution procedure is described as follows:

1. We start by looping over a subset of all insertion distances $x_c$ ranging from $0$ to $L$, in total $m$ discrete insertion distances.
2. We discretize the domain $[0,L]$ into $n$ points, meaning that $\phi_0$, $\phi_1$, $N_1$ and $\sigma(x, x_c)$ turns into an array of length $n$. 
3. Solve for the zeroth order perturbation equation $\phi_0$ using the discretized equation below

$$-\frac{D}{\Delta^2}\phi_{0, i-1}+\Big(\frac{2D}{\Delta^2}+\Sigma_a\Big)\phi_{0,i}-\frac{D}{\Delta^2}\phi_{0,i+1}=\frac{\nu}{k}\Sigma_f\phi_{0,i}$$

which can be arranged into a matrix equation of the form
$$ \mathbf{A} \mathbf{\phi_0} = \frac{\nu \Sigma_f}{k} \mathcal{1}_{n\times n} \mathbf{\phi_0} $$
This is an eigenvalue problem that can be solved for the eigenvector $\phi_0$. 

4. Use the zeroth order solution $\phi_0$ to solve for $\phi_1$ in the first order perturbation equation for some given $N_1(x)$. This is done through the following discretized equation

$$-\frac{D}{\Delta^2}\phi_{1, i-1}+\Big(\frac{2D}{\Delta^2}+\Sigma_a\Big)\phi_{1,i}-\frac{D}{\Delta^2}\phi_{1,i+1} + \textcolor{red}{N_{1,i}\sigma(x,x_c)\phi_{0,i}}=\frac{\nu}{k}\Sigma_f\phi_{1,i}$$

which once again can be arranged into a matrix equation on the form

$$\mathbf{A}\mathbf{\phi_1} + \textcolor{red}{\mathbf{S(N_1, \phi_0, x_c)}}= \frac{\nu \Sigma_f}{k} \mathcal{1}_{n\times n} \mathbf{\phi_1} $$

which can be solved for $\phi_1$ using 

$$\phi_1 = -(\mathbf{A}-\frac{\nu \Sigma_f}{k} \mathcal{1}_{n\times n})^{-1}\textcolor{red}{\mathbf{S(N_1, \phi_0, x_c)}}.$$

5. We form the final flux $\phi(x) = \phi_0(x) + \epsilon \phi_1(x)$ and the final number density $N(x) = N_0 + \epsilon N_1(x)$.

6. We use $N_1$ as a function we seek to optimize in order to minimize a cost function. This minimization will be implemented through a genetic algorithm, although it's very simple to alter the code if one wants to use standard minimization techniques such as Newton, Nelder-Mead, etc. Either way, what we are actually optimizing are the $c_\eta$'s in the Legendre expansion of $N_1$. The cost function is defined with three terms

$$ \mathrm{C}(\mathbf{c_\eta}) = \lambda_{\mathrm{PPF}} \sum_{m} (\mathrm{PPF_0} -\mathrm{PPF_1} ) +\lambda_{\mathrm{avg}} \left | \int_0^L N_1(x) dx \right | +  \lambda_{\mathrm{flux}} \left | \int_0^L \phi(x)-\phi_0(x) dx \right |.$$

The first term can be identified as the win in axial PPF between the zeroth order solution $\phi_0$ and the perturbed solution $\phi = \phi_0 + \epsilon \phi_1$ after summing over all insertions $m$. The second term punishes solutions in $N_1$ that don't average to $0$. Remember that we want to infer a design from a control rod, so we would like to use as little new material as possible when minimizing the PPF. The third contribution comes from the fact that we want to make a fair comparision between the unperturbed and perturbed system by assuring the integrated flux in the system is constant.


### Inferring a Control Rod Design

To tie everything together, we finally have to explain how we can get an idea of a control rod design given our perturbed solution $N(x)$. Although we are studying a one dimensional case, one could view it as a 3D case that's infinite in the two other spatial orthogonal axis (say $y$ and $z$), or at least that the length scale in the x direction $L_x$ fulfills $L_x << L_y, L_z$. In that case one might neglect the spatial dependencies in these directions, and it's more reasonable to infer a control rod design. We know that the volume of a cylinder is given by $\pi r^2 L$. If we let $r$ vary axially, meaning $r = r(x)$, then the number density is given by

$$ r(x) \sim \frac{N_A}{\pi r(x)^2 L} \Rightarrow r(x) \sim \sqrt{\frac{N_A}{\pi L N(x)}} $$

where $N_A$ is Avogadros constant. Hence we can get an intuition for a design by studying about how $N(x)$, or in our case since $N_0$ is constant, $N_1(x)$ varies axially. Therefore $r \sim N(x)^{-1/2}$. 

## Code Description

The code performs the optimization of a control rod design by minimizing the cost function with respect to $\mathbf{c_\eta}$. First, we initialize the code by setting relevant nuclear and spatial parameters, as well as defining relevant parameters for the Legendre expansion and the genetic algorithm. 

What follows is a construction of the matrix $A$ for the neutron diffusion equation in a multiplying medium which incorporates the effect of various insertion percentages of a control rod. 

Thereafter we have three optimization functions, firstly `invPow` that performs an inverse power iteration to find the dominant eigenvalue and its' corresponding eigenvector.`solve_phi0` and `solve_phi1` compute the zeroth-order and first-order neutron fluxes according to the defined differential equations in the background. Finally, the `objective` function defines the cost function for the genetic algorithm, and is evaluated purely from the coefficients in the Legendre polynomial expansion. 

Finally we enter the genetic algorithm that initializes a population of solutions $\mathbf{c_\nu}$'s, and selects the best individuals based on the evaluated fitness from the `objective` function. We thereafter perform crossovers and mutation in order to increase the gene pool and find new solutions. 

The fitness is then visualized dynamically during the optimization process, and each of the consistuent terms in the cost function are visualized to see their individual contributions. 

## Results

Coming soon.

## Conclusion



