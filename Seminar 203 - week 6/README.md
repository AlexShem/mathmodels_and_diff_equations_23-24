# Rod Equation: Advanced Computational Topics

## Rod Equation and Discretization

We consider the rod equation, which can be expressed in its simplified form as:

$$
\frac{\partial^2 u}{\partial t^2} - D \frac{\partial^4 u}{\partial x^2 \, \partial t^2} + C \frac{\partial^4 u}{\partial x^4} =  f(t, x),
$$

where the coefficients $D = R^2$, $C = E R^2 / \rho$, $x$ is a spatial variable, $t$ is time, $\rho>0$ is the density of the rod material, $R$ is the cross-section radius, and $E$ is Young's modulus of the material. The right part $f(t,\, x)$ is a forcing.

The discretization process introduces temporal step $\tau$ and spatial step $h$, leading to the definition of dimensionless parameters $\nu = C \tau^2/h^4$ and $\mu = D/h^2$.

## 1. Interpolation to the Final Time Moment

Usually, for the numerical simulations we set an upper integration time limit `T`. So, we expect the computational method, like, *Crank-Nicolson scheme* or a *Compact scheme* to produce the numerical solution on the time interval $t \in [0, T]$. However, by introducing the spatial step `h` and the temporal step `tau` we can only calculate the values of the unknown function $u(t, x)$ at predetermined temporal-spatial points.

For example, consider the final time of interest to be `T = 1` and the temporal step `tau = 0.3`. Then, we can only compute the function `u` at times `t = [0, 0.3, 0.6, 0.9, 1.2]`. Here, the initial conditions are set at time `t = 0` and we continue integration until the first time we exceed the desired time `T`. In this case, the total number of temporal points `Nt` needed to compute the solution is the number of termporal steps required to exceed the desired time `T` plus one for the initial time moment `t = 0`. A sample code snippet may have the following form:

```matlab
T = 1;
tau = 0.3;
Nt = ceil(T/tau) + 1;
```

![u_interpolation](rod_problems.png)

On the figure above, the solid black line represents the reference solution `u_ref`, and the blue dots correspond to the calculated values of the function `u`.

To get the value of the function `u` at the desired time `T`, we can interpolate the values using one of the available methods:

- **1-D interpolation**: Condiser each slice of the function `u` along the spatial coordinate `x` as an intependent function of time `t`, and perform `Nx` interpolations to time `T`. Combine `for` cycle and [`interp1`](https://www.mathworks.com/help/releases/R2023b/matlab/ref/interp1.html) funtction;
- **2-D interpolation**: Consider spatial and temporal variables simultaneously, and a perform gridded interpolation using [`interp2`](https://www.mathworks.com/help/releases/R2023b/matlab/ref/interp2.html).

Usually, interpolation techniques allow different `method` parameters:

- **Linear** interpolation,
- **Cubic** interpolation (based on a cubic convolution),
- **Spline** interpolation,
- and others.

Here, cubic convolution interpolation employs a kernel-based approach focusing on local smoothness, using a weighted average determined by a cubic polynomial. In contrast, cubic spline interpolation constructs a globally smooth, piecewise cubic polynomial that passes exactly through all known data points, emphasizing continuity and smoothness across the entire range, making it suitable for mathematical modeling and data fitting where precise data representation is crucial.

It is up to the user to define which method works the best in a particular problem. A code snipped may look as following:

```matlab
x = linspace(0, L, Nx + 1);

% Solve the rod equation for the given parameters and scheme
% Nx - number of spatial points;
% L - length of the rod;
% T - desired final time;
% rho, R, E - parameters of the rod;
% nu - dimensionless parameter, nu = (E*R^2/rho) * tau^2 / h^4;
% scheme - one of the computational schemes, e.g., "CN" or "555";
% u_ref, f_ref - function handles to the reference solution and the external force.
% OUTPUT:
% t - column vector of Nt time moments at which the solution is calculated;
% u - matrix of the size [Nt, Nx] with the numerical approximation of the unknown function.
[t, u] = solveRodEquation(Nx, L, T, rho, R, E, nu, scheme, u_ref, f_ref);

u_numerical = interp2(x, t, u, x, T, "spline");
```

