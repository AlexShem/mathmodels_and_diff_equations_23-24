# Practical Session: Numerical Solutions to the Rod Equation

## Overview

In this session, we will investigate the numerical solutions of the rod equation using different finite difference schemes. The primary focus will be on the standard Crank-Nicolson scheme and the development of compact schemes, aiming to understand their implementation and effectiveness in solving the rod equation numerically.

## Goals

1. Understand the discretization of the rod equation and its implications.
2. Implement the Crank-Nicolson scheme and a compact scheme in MATLAB.
3. Analyze the accuracy and stability of the numerical solutions.

## Section 1: Theoretical Background

### Rod Equation and Discretization

We consider the rod equation, which can be expressed in its simplified form as:

$$
\frac{\partial^2 u}{\partial t^2} - D \frac{\partial^4 u}{\partial x^2 \, \partial t^2} + C \frac{\partial^4 u}{\partial x^4} =  f(t, x),
$$
where the coefficients $D = R^2$, $C = E R^2 / \rho$, $x$ is a spatial variable, $t$ is time, $\rho>0$ is the density of the rod material, $R$ is the cross-section radius, and $E$ is Young's modulus of the material. The right part $f(t,\, x)$ is a forcing.

The discretization process introduces temporal step $\tau$ and spatial step $h$, leading to the definition of dimensionless parameters $\nu = C \tau^2/h^4$ and $\mu = D/h^2$.

### Boundary Conditions

We apply periodic boundary conditions, implying the equivalence of values on the left and right borders of the spatial domain.

## Section 2: Numerical Problem Statement

### Grid and Time Step Determination

Given a rod of length $L$, we compute the spatial step $h = L/N$ for $N$ spatial points. The temporal step $\tau$ is then calculated based on $h$ and the fixed dimensionless parameter $\nu$.

### Reference Solutions

We will consider two reference solutions $u^\mathrm{ref}$ for our numerical experiments, with the aim to compare them to the obtained numerical solutions $u$ and evaluate the error norms.

## Section 3: MATLAB Implementation

### Task 3.1: Discrete Grid Setup

Establish the discrete grid in MATLAB, considering the number of spatial points $N$ and temporal steps $M$.

### Task 3.2: Coefficient Calculation

Calculate the coefficients for the CN scheme and the compact scheme as per the given formulas.

**Code Snippet**:

```matlab
% Coefficients for CN Scheme
alpha = 1 + 3*nu + 2*mu;
% Left part
a = -(2*nu + mu) / alpha;
b = -(2 + 4*mu) / alpha;
c = 2*mu / alpha;
d = 0;
e = nu/2 / alpha;
% Right part
p0 = 0;
p1 = 0;
q0 = 0;
q1 = 0;
r0 = 0;
r1 = tau^2/2 / alpha;
```

```matlab
% Coefficients for "5-5-5" Scheme
alpha = 72*(41 + 30*nu + 90*mu);
% Left part
a = 96*(7 - 15*nu - 30*mu) / alpha;
b = -144*(41 - 150*nu + 90*mu) / alpha;
c = -192*(7 + 75*nu - 30*mu) / alpha;
d = -24*(1 - 150*nu - 30*mu) / alpha;
e = 12*(1 + 30*nu - 30*mu) / alpha;
% Right part
p0 = 10*tau^2 / alpha;
q0 = 560*tau^2 / alpha;
r0 = 2460*tau^2 / alpha;
p1 = tau^2 / alpha;
q1 = 56*tau^2 / alpha;
r1 = 246*tau^2 / alpha;
```

## Section 4: Numerical Schemes Comparison

### Task 4.1: Define the Problem Parameters

Consider a scenario with the steel rod with the following physical parameters:

```matlab
rho = 7900;
R = 10^-3;
E = 210e9;
L = 2*pi;
```

Set the integaration time to `T = 1` and the number of spatial points to be `Nx = 50`. Also, fix the dimensionless parameter $\nu = C \tau^2/h^4$ to be equal to $0.1$. Calculate the corresponding time step $\tau$.

```matlab
T = 1;
Nx = 100;
nu = 0.1;

x = linspace(0, L, Nx + 1);
h = x(2) - x(1);
D = R^2;
C = E*R^2/rho;
mu = D/h^2;

% Secondary parameters
tau = sqrt(nu*h^4 / C);
Nt = ceil(T/tau) + 1; % Including the t = 0 and the first step beyond T
T_final = tau * (Nt-1);
t = linspace(0, T_final, Nt);
```

Finally, define the reference solution and the corresponding forcing:

```matlab
% Reference soltuion
u_ref = @(t, x) sin(x + t*sqrt(C/(D+1)));
f_ref = @(t, x) zeros(size(x));
```

### Task 4.2: Implementing the CN and "5-5-5" Schemes

Implement the CN and "5-5-5" schemes in MATLAB and solve the rod equation numerically. Start by defining the transition matrices for the left- and right-hand sides.

**Code Snippet**:

```matlab
U_next = eye(Nx) + diag(a*ones(Nx-1, 1), 1) + diag(a*ones(Nx-1, 1), -1) + ...
    diag(e*ones(Nx-2, 1), 2) + diag(e*ones(Nx-2, 1), -2);
U_now = b*eye(Nx) + diag(c*ones(Nx-1, 1), 1) + diag(c*ones(Nx-1, 1), -1) + ...
    diag(d*ones(Nx-2, 1), 2) + diag(d*ones(Nx-2, 1), -2);

F_next = r1*eye(Nx) + diag(q1*ones(Nx-1, 1), 1) + diag(q1*ones(Nx-1, 1), -1) + ...
    diag(p1*ones(Nx-2, 1), 2) + diag(p1*ones(Nx-2, 1), -2);
F_now = r0*eye(Nx) + diag(q0*ones(Nx-1, 1), 1) + diag(q0*ones(Nx-1, 1), -1) + ...
    diag(p0*ones(Nx-2, 1), 2) + diag(p0*ones(Nx-2, 1), -2);

% Periodic boundary conditions
U_next(1, end-1:end) = [e a];
U_next(2, end) = e;
U_next(end-1, 1) = e;
U_next(end, 1:2) = [a e];

U_now(1, end-1:end) = [d c];
U_now(2, end) = d;
U_now(end-1, 1) = d;
U_now(end, 1:2) = [c d];

F_next(1, end-1:end) = [p1 q1];
F_next(2, end) = p1;
F_next(end-1, 1) = p1;
F_next(end, 1:2) = [q1 p1];

F_now(1, end-1:end) = [p0 q0];
F_now(2, end) = p0;
F_now(end-1, 1) = p0;
F_now(end, 1:2) = [q0 p0];

U_prev = U_next;
F_prev = F_next;
```

Now, implement the numerical method.

```matlab
% Integration
u = zeros(Nt, Nx);
u(1, :) = u_ref(0, x(1:end-1));
u(2, :) = u_ref(tau, x(1:end-1));
[x_gr, t_gr] = meshgrid(x(1:end-1), t);

f = f_ref(t_gr, x_gr);

for k = 3 : Nt
    rhs = -U_now*u(k-1, :).' - U_prev*u(k-2, :).' + ...
        F_next*f(k, :).' + F_now*f(k-1, :).' + F_prev*f(k-2, :).';
    u(k, :) = U_next \ rhs;
end
u = [u, u(:, 1)];
```

## Section 5: Error Analysis

### Task 5.1: Error Norm Calculation

Calculate the $C$ and $L^2$ norms of the error between the reference solution and the numerical solution obtained from each scheme at the final time moment. To do so, run the numerical simulations for a selected scheme for several values of $N$ by running `rod_integration` function (created from the code above). Each time, calculate and store the errors.

```matlab
scheme = "CN";
Nx_set = [12, 25, 50, 100];

C_norm = zeros(length(Nx_set), 1);
L2_norm = zeros(length(Nx_set), 1);

for n = 1:length(Nx_set)
    Nx = Nx_set(n);
    [u, T_final] = rod_integration(Nx, scheme, ...) % and other parameters, if necessary

    diff_u = u_ref(T_final, x) - u(end, :)
    C_norm(n) = max(abs(diff_u));
    L2_norm(n) = sqrt(trapz(x, diff_u.^2));
end
```

### Task 5.2: Visualization

Visualize the error norms and the numerical solutions compared to the reference solutions.

**Code Snippet**:

```matlab
figure(1);
loglog(Nx_set, C_norm);
xlabel('N_x');
title('C norm');

figure(2);
loglog(Nx_set, L2_norm);
xlabel('N_x');
title('L^2 norm');
```

## Conclusion

Summarize the findings from the practical session, highlighting the differences in accuracy between the schemes.
