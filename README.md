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
