# <img src="mmde_logo.png" alt="mmde_logo" width="35" height="35"> Mathematical Models and Differential Equations

# 📅 Schedule 

Sometimes on Thursdays
- **Standard time**: 13:00 - 14:20
    - Zoom link: [https://unil.zoom.us/j/96715266118](https://unil.zoom.us/j/99765979512)

# 📢 Announcements

* MATLAB class on September 21th.

# 📚 Content

## Seminar 1: _Transport and Continuity Equations_

September 21st, 2023

- Characteristics of the equation: $d_t x(t) = u(t, x)$
- Solution when $u(t, x) = ax + b$
- Development of formulas
- Change of the density $\rho(t, x(t))$ along characteristics $x(t)$
- Numeric calculations using `ode45()`
- Visualizations of the results

## Seminar 2: _Fourier Method for Diffusion Equation_

September 28th, 2023

$$\partial_t u(t, x) = D \partial_x^2 u(t, x)$$

- Intution behind the equation contruction
- Dirichlet boundary conditions
- Second order spatial derivative
    - Eigenvalues and eigenfunctions
    - Border conditions satisfaction
- Development of the formulas

## Seminar 3: _Fourier Method for Diffusion Equation (continued)_

October 5th, 2023

- Solution in the special form
- Development of the final solution form
- Scalar product in the $L_2$ space
- Code of the Fourier method
    - Algorithm
    - Visualisation (with animation)

## Seminar 4: _Explicit Euler Scheme for Diffusion Equation_

October 12th, 2023

- Numerical approximation of the first and second derivatives
- Statement of the explicit Euler scheme for the diffusion equation
- Graphical form of the scheme
- Matrix form of the scheme
- Dirichlet border conditions in the matrix form

## Seminar 5: _Euler Schemes and Stability Analysis for Diffusion Equation in Matlab_

November 9th, 2023

- Implementation of the explicit Euler method in Matlab.
- First time step, when it crossed the boundary `T`.
- Filling in the three-diagonal matrix of the system using `eye()` and `diag()`.
- Integrating the system over time using `for` loop.
- Visualization of the solution as a surface over the space-time plane `surf(x, y, u)`.
- Implementation of Neumann boundary conditions $\partial_x u(t, 0) = 0$ for the explicit Euler scheme.
- Animation of the solution over time using `drawnow`.
- Implementation of the implicit Euler scheme.
- Experimental comparison of the stability regions of both schemes depending on $\nu = D \tau / h^2$.

## Seminar 6: _Constructing and Implementing Compact Schemes for the Heat Equation_

November 23th, 2023

- Development of a compact scheme for the diffusion equation.
- Designing a scheme template and normalization of its coefficients.
- Selection and application of test functions into the compact scheme equation to derive the coefficient equations.
- Symbolic solution of the coefficient system using `solve()`.
- Computing the heat equation (diffusion) solution using the derived compact scheme.
- Implementation of periodic boundary conditions in the scheme.

## Seminar 7: _Enhancing Compact Schemes for the Diffusion Equation with External Force_

December 7th, 2023

- **Introduction to External Force Integration in Diffusion Equations**:
    - Modifying existing numerical schemes to include external force terms.
- **Defining and Implementing the External Force Function**:
    - Creation of a MATLAB function to define the external force $f(t, x)$.
    - Integration of this force function into the numerical scheme for the diffusion equation.
- **Establishing a True Solution for Comparison**:
    - Definition of a true solution $u_{\text{true}}(t, x) = \sin(t) \cdot \cos(x)$ for the diffusion equation with an external force.
    - MATLAB implementation of the true solution for accuracy comparison.
- **Derivation and Implementation of the Corresponding Force**:
    - Analytical derivation of the force function corresponding to the chosen true solution.
    - Implementing this derived force function in MATLAB.
- **Numerical and True Solution Comparison**:
    - Visualization and comparison of the numerical solution with the true solution at various time steps.
    - MATLAB code for plotting and analyzing the solution differences.
- **Grid Size Analysis and Accuracy Assessment**:
    - Evaluation of the numerical solution accuracy for different grid sizes.
    - Keeping the Courant parameter $\nu$ constant during this analysis.
    - Plotting and interpreting the $C$ and $L^2$​ norms of the difference between numerical and true solutions.
- **Hands-On MATLAB Coding Session**:
    - Analysis of solution accuracy against the true solution for various grid sizes.

## Seminar 201: _Numerical Solutions to the Rod Equation_

January 18th, 2024

- Investigation of numerical solutions for the rod equation using Crank-Nicolson and compact schemes.
- Implementation and analysis of the Crank-Nicolson and "5-5-5" compact schemes in MATLAB to solve the rod equation.
- Periodic boundary conditions for the rod equation.
- Comparison of reference solutions with the obtained numerical solutions to evaluate accuracy and stability.
- Error analysis through calculating error norms at the final time moment for different spatial points using the implemented numerical schemes.

## Seminar 202: _Modular Implementation and Analysis of Numerical Schemes for the Rod Equation in Matlab_

February 1st, 2024

- Creation of a modular Matlab code for simulating the rod equation using Crank-Nicolson (CN) and "5-5-5" schemes.
- Explanation of the computational process involving spatial and temporal discretization parameters.
- Application and comparison of reference solutions with numerical solutions for accuracy assessment.
- Analysis of errors using C and L^2 norms for different numbers of spatial steps.
- Visualization and evaluation of the numerical scheme's performance through graphs plotting error norms against spatial discretization points.
