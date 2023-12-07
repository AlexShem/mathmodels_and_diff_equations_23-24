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
