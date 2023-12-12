# Exercise: Enhancing a Compact Scheme for the Diffusion Equation with External Force

# Overview

This exercise involves modifying a MATLAB script that implements a compact scheme for solving the diffusion equation. The task is to incorporate an external force, compare the numerical solution with a true solution, and analyze the accuracy for various grid sizes.

## Goals

1. **Incorporate External Force**: Enhance the MATLAB code to account for an external force $f(t, x)$.
2. **True Solution and Force Determination**: Establish a true (reference) solution $u_{\text{true}}(t, x) = \sin(t) \cdot \cos(x)$ and deduce the corresponding force $f(t, x)$.
3. **Solution Comparison**: Compare the solution from the numerical method with the true solution.
4. **Grid Size Analysis**: Analyze the solutions for multiple grid sizes while keeping the Courant parameter $\nu = \frac{D \tau}{h}$ fixed. Plot the $C$ and $L_2$ norms of the difference between the obtained and true solutions.

## Section 1: Incorporating External Force

Modify the provided MATLAB code to include an external force term in the diffusion equation. The force function $f(t, x)$ should be defined and integrated into the numerical scheme.

### Task 1.1: Define the External Force

Define a new function for the external force $f(t, x)$. This function should be able to accept time $t$ and spatial coordinate $x$ as inputs.

**Code Snippet**:

```matlab
function f_val = external_force(t, x)
    % Define the external force function here
    f_val = ...; % Replace with the actual function
end
```

### Task 1.2: Update the Numerical Scheme

Integrate the external force into the numerical scheme. Ensure the force term is appropriately added to the scheme's update step.

**Code Modification**:

```matlab
% Inside the integration loop
for k = 2 : Nt
    % Calculate the force term for the current time step
    % Note: you might need to do that for several time steps
    % depending on the number of temporal layers of the scheme
    f_term = external_force((k-1)*tau, x(1:end-1));
    % Update U(k, :) to include the force term
    U(k, :) = -U_next \ (U_now * U(k - 1, :).' + f_term);
end
```

## Section 2: True Solution and Force Determination

### Task 2.1: Define the True Solution

Define the true solution $u_{\text{true}}(t, x) = \sin(t) \cdot \cos(x)$.

**Code Snippet**:

```matlab
function u_true_val = u_true(t, x)
    u_true_val = sin(t) .* cos(x);
end
```

### Task 2.2: Determine the Corresponding Force

Derive the formula for the force $f(t, x)$ that corresponds to the given true solution. Provide the derivation in the report and implement the function in MATLAB.

## Section 3: Solution Comparison

Compare the numerical solution with the true solution at various time steps. Plot the solutions for visual comparison.

### Task 3.1: Compare at Specific Time Steps

**Code Snippet**:

```matlab
% Comparison loop
for k = [selected time steps]
    % Calculate true solution
    u_true_k = u_true((k-1)*tau, x);
    % Plot numerical and true solutions
    figure;
    plot(x, U(k, :), 'r', x, u_true_k, 'b--');
    legend('Numerical Solution', 'True Solution');
    title(['t = ', num2str((k-1)*tau)]);
end
```

## Section 4: Grid Size Analysis

Analyze the accuracy of the numerical solution for different grid sizes. Keep the Courant parameter $\nu$ fixed and plot the $C$ and $L_2$ norms of the differences.

### Task 4.1: Perform Analysis for Various Grid Sizes

**Code Snippet**:

```matlab
% Define different grid sizes
grid_sizes = [50, 100, 200, 400, ...];

% Loop over grid sizes
for Nx = grid_sizes
    % Update other parameters based on Nx
    % ...
    % Perform the simulation
    % ...
    % Calculate and plot C and L2 norms
    % ...
end
```

### Task 4.2: Plotting Norms

Plot the $C$ and $L_2$ norms for each grid size. Use logarithmic scaling if necessary to better visualize the results.

## Coefficients of Computational Schemes

This section outlines the coefficients used in three computational schemes: the Explicit Euler Scheme, the Crank-Nicolson Scheme, and the Compact Scheme. These coefficients are crucial for the stability and accuracy of numerical solutions to differential equations.

### 1. Explicit Euler Scheme (Euler\_ex)

The Explicit Euler Scheme is a straightforward time integration method used for solving ordinary differential equations. It is a first-order method.

#### Coefficients:

*   Left Part:
    *   $a = 0$
    *   $b = -1 + 2\nu$
    *   $c = -\nu$
*   Right Part:
    *   $p1 = 0$
    *   $p2 = 0$
    *   $q1 = \tau$
    *   $q2 = 0$

#### MATLAB Code Snippet:

```matlab
% Explicit Euler Scheme Coefficients
a = 0;
b = -1 + 2*nu;
c = -nu;
p1 = 0;
p2 = 0;
q1 = tau;
q2 = 0;
```

### 2. Crank-Nicolson Scheme (cn)

The Crank-Nicolson Scheme is a popular second-order method for numerical solution of partial differential equations.

#### Coefficients:

*   Left Part:
    *   $a = -\frac{\nu}{2\alpha}$
    *   $b = \frac{-1 + \nu}{\alpha}$
    *   $c = -\frac{\nu}{2\alpha}$
*   Right Part:
    *   $p1 = 0$
    *   $p2 = 0$
    *   $q1 = \frac{\tau}{2\alpha}$
    *   $q2 = \frac{\tau}{2\alpha}$
*   Where $\alpha = 1 + \nu$

#### MATLAB Code Snippet:

```matlab
% Crank-Nicolson Scheme Coefficients
alpha = 1 + nu;
a = -nu/2/alpha;
b = (-1 + nu)/alpha;
c = -nu/2/alpha;
p1 = 0;
p2 = 0;
q1 = tau/2/alpha;
q2 = tau/2/alpha;
```

### 3. Compact Scheme

The Compact Scheme is designed for high accuracy and efficiency in solving differential equations.

#### Coefficients:

*   Left Part:
    *   $a = \frac{2(1 - 6\nu)}{\alpha}$
    *   $b = \frac{-4(5 - 6\nu)}{\alpha}$
    *   $c = \frac{-2(1 + 6\nu)}{\alpha}$
*   Right Part:
    *   $p1 = \frac{\tau}{\alpha}$​
    *   $p2 = \frac{\tau}{\alpha}$​
    *   $q1 = \frac{10\tau}{\alpha}$
    *   $q2 = \frac{10\tau}{\alpha}$
*   Where $\alpha = 4(6\nu + 5)$

#### MATLAB Code Snippet:

```matlab
% Compact Scheme Coefficients
alpha = 4*(6*nu + 5);
a = 2*(1 - 6*nu)/alpha;
b = -4*(5 - 6*nu)/alpha;
c = -2*(1 + 6*nu)/alpha;
p1 = tau/alpha;
p2 = tau/alpha;
q1 = 10*tau/alpha;
q2 = 10*tau/alpha;
```
