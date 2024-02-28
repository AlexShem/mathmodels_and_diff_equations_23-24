%% Primary Parameters Definition
% L - Length of the rod in meters.
% T - Desired final time in seconds.
% rho - Density of the rod material in kg/m^3.
% R - Radius of the rod in meters.
% E - Young's modulus of the rod material in Pascals.
% nu - Dimensionless parameter, calculated as (E*R^2/rho) * tau^2 / h^4, governing the equation's behavior.
% scheme - Specifies the computational scheme used for the simulation, e.g., "CN" (Crank-Nicolson) or "555".
% omega - Time frequency in radians per second, related to the forcing function.
% Nx - Number of spatial points used in the discretization of the rod.
L = pi;
T = 1;
rho = 7900;
R = 1e-2;
E = 210e9;
nu = 0.1;
omega = pi;
Nx = 45;

%% Selection of Scheme and Boundary Conditions
% scheme - The numerical scheme used for solving the differential equation. Choices include "CN" and "555".
% boundary_conditions - Specifies the type of boundary conditions applied, e.g., "u0u1" or "u0u2".
% approximation_type - Defines the approximation approach used, options are "naive" or "compact".
scheme = "555";
boundary_conditions = "u0u2";
approximation_type = "compact";

%% Calculation of Secondary Parameters
% D, C - Derived parameters based on the physical properties of the rod.
% x - Spatial grid points along the length of the rod.
% h - Spatial step size calculated from the grid points.
% mu - A derived parameter for the model, related to the discretization and material properties.
% tau - Temporal step size calculated from the nu parameter and physical constants.
% Nt - Number of temporal points calculated based on the final time and temporal step size.
% t - Temporal grid points.
D = R^2;
C = E*R^2/rho;
x = linspace(0, L, Nx);
h = x(2) - x(1);
mu = D / h^2;
tau = sqrt(nu*h^4/C);
Nt = ceil(T/tau) + 1;
t = (0 : tau : tau*(Nt-1)).';

%% Reference Solutions and External Force Functions
% u_ref - Reference solution function handle, dependent on time and space.
% f_ref - External force function handle, dependent on time and space.
% These functions are defined based on the selected boundary conditions.
if boundary_conditions == "u0u1"
    c_t = @(t) cos(omega * t);
    u_ref = @(t, x) c_t(t) .* sin(x).^2;
    f_ref = @(t, x) cos(omega.*t).*((omega.^2.*cos(x.*2.0))./2.0-omega.^2./2.0-C.*cos(x.*2.0).*8.0+D.*omega.^2.*cos(x.*2.0).*2.0);
elseif boundary_conditions == "u0u2"
    c_t = @(t) cos(omega * t);
    u_ref = @(t, x) c_t(t) .* sin(x);
    f_ref = @(t, x) -cos(omega.*t).*sin(x).*(-C+D.*omega.^2+omega.^2);
end

%% Coefficients of the Computational Scheme
% Retrieves the coefficients required for the selected numerical scheme and boundary conditions.
[~, a, b, c, d, e, p0, p1, q0, q1, r0, r1] = getCoefficients(scheme, nu, D, h, tau);

%% Construction of Computational Matrices
% This function constructs matrices based on the Nx, nu, mu, C, D, h, tau, boundary_conditions, approximation_type, and scheme coefficients.
% These matrices are used in the numerical solution of the differential equation.
[U_next, U_now, U_prev, F_next, F_now, F_prev] = getMatrices(Nx, nu, mu, C, D, h, tau, boundary_conditions, approximation_type, a, b, c, d, e, p0, p1, q0, q1, r0, r1);

%% Initialization of Solution and External Force Arrays
% u - Array to store the solution at each spatial and temporal point.
% u(1, :) and u(2, :) - Initial conditions applied at t=0 and t=tau, respectively.
% u_star - Reference solution array for comparison.
% f - External force array applied throughout the simulation.
u = zeros(Nt, Nx); % Initialize the solution matrix with zeros.
u(1, :) = u_ref(0, x); % Set initial condition at t=0.
u(2, :) = u_ref(tau, x); % Set initial condition at t=tau for stability.

[x_gr, t_gr] = meshgrid(x, t); % Create mesh grids for space and time.
u_star = u_ref(t_gr, x_gr); % Calculate the reference solution over the grid.
f = f_ref(t_gr, x_gr); % Calculate the external force over the grid.

%% Time-Stepping Loop for Numerical Integration
% Iteratively solves the differential equation using the computed matrices and initial conditions.
for k = 3:Nt
    % Calculate the right-hand side (rhs) of the equation based on the previous two time steps and external force.
    rhs = -U_now*u(k-1, :).' - U_prev*u(k-2, :).' + F_next*f(k, :).' + F_now*f(k-1, :).' + F_prev*f(k-2, :).';
    % Solve for the current time step using the rhs and U_next matrix.
    u(k, :) = U_next \ rhs;
end

%% Visualization of Solution via Animation
% Creates a figure to visually compare the numerical solution (u) and the reference solution (u_star) over time.
figure(1)
for k = [1 : 100 : Nt, Nt] % Iterates over selected time steps for visualization.
    plot(x, u(k, :), '-b'); % Plot numerical solution.
    hold on;
    plot(x, u_star(k, :), '--r'); % Plot reference solution.
    hold off;
    xlabel('$x$', 'FontSize', 16, 'Interpreter', "latex"); % Label for the x-axis.
    title(['$t = $', num2str(round(t(k), 3))], 'FontSize', 16, 'Interpreter', "latex"); % Dynamic title showing current time.
    axis([0, L, -1.1, 1.1]); % Set axis limits.
    drawnow limitrate nocallbacks; % Update the figure.
end

%% Calculation and Visualization of Error
% Computes the maximum error at each time step and visualizes it in a semi-log plot.
C_norm_t = max(abs(u - u_star), [], 2); % Calculate the maximum error over all spatial points at each time step.
figure(3); % Create a new figure for the error plot.
semilogy(t, C_norm_t); % Plot the error over time in a semi-log plot.
hold on;
