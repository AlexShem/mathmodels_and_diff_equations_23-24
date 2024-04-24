# Optimization Techniques

## Introduction

In this session, we will explore two optimization techniques - **Gradient Descent** and **Particle Swarm Optimization (PSO)** - to find the minimum of a 2-dimensional function in a specified area. The objective is to gain a deeper understanding of these methods, implement them, and discuss their benefits and limitations.

### Task Description

Your task is to minimize a given 2-dimensional function `f(x, y)` and find its global minimum within the specified area. You will implement this using two methods: Gradient Descent and Particle Swarm Optimization. After implementing these methods, you should compare their performance and discuss the benefits and limitations of each approach.

### Function to Minimize

Consider the function defined as:

$$
f(x, y) = \exp\left(-\frac{x^2 + y^2}{10}\right) \cdot \sin(x) \cdot \cos(y)
$$

This function has multiple local minima and maxima, making it an interesting case for optimization.

### Area of Interest

The search will be conducted in the area defined by:

- Lower bounds: $[-4, -4]$
- Upper bounds: $[4, 4]$

## Methods to Use

### 1. Gradient Descent

Gradient Descent is a first-order iterative optimization algorithm for finding the minimum of a function. Here is a snippet to guide you in implementing Gradient Descent for our function:

```matlab
%% Gradient descent: Initialization
x0 = [1.5; 1.5];  % Initial point
[xmin_g, fmin_g, niter_g, path_g] = grad_descent(x0, fx, 1e-8, 1000);

%% Gradient descent: Visualization
figure(1)
hold on;
plot3(path_g(1,:), path_g(2,:), f(path_g(1,:), path_g(2,:)), ...
    '-*r', 'LineWidth', 2, 'MarkerSize', 12)
hold off;
```

### 2. Particle Swarm Optimization (PSO)

PSO is a computational method that optimizes a problem by iteratively improving a candidate solution with regard to a given measure of quality. Implement the given PSO function and observe the optimization process:

```matlab
%% Particle Swarm Optimization: Visualization and Implementation
figure(3)
fsurf(f, [lb(1) ub(1) lb(2) ub(2)], 'ShowContours', 'on', 'FaceAlpha', .3)

[xmin_s, fmin_s] = particle_swarm(fx, f, lb, ub, 20, 50);
```

## Analysis

After implementing both methods, analyze the results:

- Which method converged faster?
- How did the choice of initial points affect the results?
- Compare the robustness of each method against local minima.

## Assignment Submission

Please submit a report that includes:

- Your implementation code for both methods.
- A discussion on the behavior of both optimization methods, including their convergence behaviors.
- Screenshots of the plots generated during the optimization processes.
- A comparative analysis of the benefits and limitations of Gradient Descent and Particle Swarm Optimization based on your observations.

Feel free to use the provided code snippets to fill in during your practical session. Explore different initial conditions and parameters to fully understand the behavior of these optimization techniques.

---

## Implementation of Gradient Descent

### Overview of Gradient Descent

Gradient Descent is a popular optimization algorithm used in many fields, particularly in machine learning and optimization of complex mathematical functions. It's a first-order iterative optimization algorithm for finding the minimum of a function. The idea is straightforward: take repetitive steps in the opposite direction of the gradient (rise) of the function at the current point because this is the direction of steepest descent.

### Approach to Implementing Gradient Descent

The basic steps in a Gradient Descent algorithm include:

1. **Initialization**: Start with an initial guess `x0` for the location of the minimum.
2. **Gradient Calculation**: Compute the gradient (derivative) of the function at this point.
3. **Update Step**: Move in the direction opposite to the gradient by a step size (also known as learning rate) to get to the next point.
4. **Convergence Check**: Repeat the process until the change in the function value is less than a specified tolerance or a maximum number of iterations is reached.

These steps make Gradient Descent particularly suitable for problems where the objective function is differentiable.

### Detailed Implementation Steps

Here's how you can implement the Gradient Descent method based on the provided framework:

1. **Function Signature and Parameters**: Set up the function to take initial parameters like initial guess `x0`, the function to minimize `f`, tolerance `tol`, and maximum number of iterations `maxiter`.

2. **Initial Setup**: Initialize variables such as the step path and iteration counter. Prepare to store each intermediate state to visualize the path later.

3. **Main Iteration Loop**:

    - **Gradient Calculation**: Use central differences to compute an approximation of the gradient at the current point.
    - **Step Update**: Update the position by taking a step in the opposite direction of the gradient. Adjust the step size if the new function value does not improve.
    - **Convergence Check**: Exit the loop if the change in the function value is less than the tolerance or the maximum iterations are reached.

4. **Post-Iteration**: Handle any post-iteration tasks like warnings if convergence was not achieved.

### MATLAB Code for Gradient Descent

Below is a snippet you can use to fill in during your practical session:

```matlab
% Initialization and setup
x0 = [1.5; 1.5];  % Adjust the initial guess as needed
tol = 1e-8;       % Tolerance for convergence
maxiter = 1000;   % Maximum number of iterations

% Call the gradient descent function
[xmin_g, fmin_g, niter_g, path_g] = grad_descent(x0, fx, tol, maxiter);

% Visualization of the path
figure(1)
hold on;
plot3(path_g(1,:), path_g(2,:), arrayfun(@(x, y) f(x, y), path_g(1,:), path_g(2,:)), ...
    '-*r', 'LineWidth', 2, 'MarkerSize', 12)
hold off;
```

### Gradient Calculation Function `grad_f`

This is a critical component. It calculates the gradient using central differences, which provides a better approximation than forward differences:

```matlab
function grad = grad_f(x, f)
d = length(x);
h = 1e-8;
grad = arrayfun(@(ind) (f(x + h*((1:d).' == ind)) - f(x - h*((1:d).' == ind)))/(2*h), (1:d).');
grad = grad / norm(grad, 2);  % Normalize to avoid too large steps
end
```

By following this approach and utilizing the provided code snippets, you can implement an effective gradient descent optimizer for the function given in your assignment.

---

## Implementation of Particle Swarm Optimization (PSO)

### Overview of Particle Swarm Optimization

Particle Swarm Optimization (PSO) is a robust stochastic optimization technique based on the movement and intelligence of swarms. PSO applies the concept of social interaction to solve optimization problems in a manner that mimics the behavior of biological populations.

### Approach to Implementing PSO

PSO initializes a group of random particles (solutions) and then searches for optima by updating generations. Each particle adjusts its trajectory towards its personal best location and the group’s best location according to simple mathematical formulae over the particle's position and velocity. Each movement of a particle is influenced by its local best known position and is also guided toward the best known positions in the search-space, which are updated as better positions are found by other particles.

### PSO: Detailed Implementation Steps

Here is how you can implement Particle Swarm Optimization for minimizing our given function:

1. **Initialization**:

   - Initialize the swarm with random positions and velocities within the defined bounds.
   - Evaluate the fitness of each particle at its initial position.

2. **Iteration Loop**:

   - Update the velocity of each particle based on both its own best known position and the swarm’s best known position.
   - Update the position of each particle according to its new velocity and ensure it remains within bounds.
   - Re-evaluate the fitness of each particle at its new position.
   - Update the personal best and global best positions based on the new fitness.

3. **Termination**:
   - The algorithm terminates when a stopping criterion is met, typically a maximum number of iterations or a sufficient solution fitness.

### MATLAB Code for Particle Swarm Optimization

Below is a detailed explanation and corresponding MATLAB code snippets for implementing PSO:

```matlab
% Particle Swarm Optimization: Initialization
S = 20;  % Number of particles
niter = 50;  % Number of iterations
lb = [-4, -4];  % Lower bounds
ub = [4, 4];  % Upper bounds

% Call the particle swarm function
[xmin_pso, fmin_pso] = particle_swarm(fx, f, lb, ub, S, niter);

% Visualization of the swarm's progress
figure(3)
fsurf(f, [lb(1) ub(1) lb(2) ub(2)], 'ShowContours', 'on', 'FaceAlpha', .3)
xlabel('x');
ylabel('y');
```

Below is a MATLAB snippet that implements the pseudo-code for Particle Swarm Optimization (PSO) described above. This snippet provides a framework that can be filled in with specific function and parameter details as required for your optimization problem.

```matlab
% Define the function to minimize, domain boundaries and parameters
f = @(x) yourObjectiveFunction(x);
lb = [-4, -4]; % Lower bounds of the search space
ub = [4, 4];   % Upper bounds of the search space
S = 20;        % Number of particles
niter = 100;   % Maximum number of iterations

% PSO parameters
omega = 0.4;   % Inertia weight
phip = 1.2;    % Cognitive coefficient
phig = 0.9;    % Social coefficient

% Initialize the particles
n = length(lb); % Number of dimensions
X = rand(S, n);
% Scale and shift to the propper domain
for i = 1 : n
    X(:, i) = X(:, i) * (ub(i) - lb(i)) + lb(i);
end
P = X;          % Best known positions of particles
% Initial velocity
V = rand(S, n);
for i = 1 : n
    span = ub(i) - lb(i);
    V(:, i) = span*2*(V(:, i) -.5);
end

% Evaluate the fitness at the initial positions
Fp = arrayfun(@(i) f(P(i,:)), 1:S);
[Fg, g_ind] = min(Fp);
G = P(g_ind, :);  % Global best position

% Main PSO loop
iter = 0;
while iter < niter
    for i = 1:S
        for d = 1:n
            rp = rand;
            rg = rand;
            % Update velocity
            V(i,d) = omega * V(i,d) + phip * rp * (P(i,d) - X(i,d)) + phig * rg * (G(d) - X(i,d));
        end
        % Update position
        X(i,:) = X(i,:) + V(i,:);
        % Ensure particles stay within bounds
        X(i,:) = max(X(i,:), lb);
        X(i,:) = min(X(i,:), ub);

        % Evaluate the new fitness
        new_fitness = f(X(i,:));
        if new_fitness < Fp(i)
            P(i,:) = X(i,:);
            Fp(i) = new_fitness;
            % Update the global best if necessary
            if new_fitness < Fg
                G = X(i,:);
                Fg = new_fitness;
            end
        end
    end
    iter = iter + 1;
end

% Output the results
xmin = G;
fmin = Fg;
```

### PSO: Explanation

1. **Initialization**: Particles' positions and velocities are initialized within specified bounds. Each particle's best known position is set to its initial position. The global best known position is initialized based on the best initial fitness.

2. **Velocity and Position Update**: For each particle and each dimension, velocities are updated using the inertia weight, cognitive component (based on the particle's best position), and social component (based on the global best position). Positions are then updated based on the new velocities and clamped within the search space boundaries.

3. **Fitness Evaluation and Updating Best Positions**: The fitness of each particle at its new position is evaluated. If this new fitness is better than its previous best, the particle's best known position is updated. If a particle's new best position is better than the global best, the global best is also updated.

4. **Termination**: The loop continues until a termination criterion (in this case, the maximum number of iterations) is met.

This MATLAB code snippet effectively translates the given pseudo-code into a working implementation of PSO that can be adapted for specific optimization problems by defining the objective function and adjusting parameters as needed.
