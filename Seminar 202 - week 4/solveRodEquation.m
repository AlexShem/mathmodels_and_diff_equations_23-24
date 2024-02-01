function [t, u] = solveRodEquation(Nx, L, T, rho, R, E, nu, scheme, u_ref, f_ref)
%% Calculate primary parameters
D = R^2;
C = E*R^2/rho;
x = linspace(0, L, Nx + 1);
h = x(2) - x(1); % Spatial step size
tau = sqrt(nu*h^4 / C); % Time step size, derived from nu
Nt = ceil(T/tau) + 1; % Number of time steps
t = 0 : tau : tau*(Nt-1);

%% Initial conditions and External force
u = zeros(Nt, Nx); % Initialize solution matrix
u(1, :) = u_ref(0, x(1:end-1)); % Apply initial condition at t=0
u(2, :) = u_ref(tau, x(1:end-1)); % Apply initial condition at t=tau

[x_gr, t_gr] = meshgrid(x(1:end-1), t);
f = f_ref(t_gr, x_gr);

%% Coefficients for the computational scheme
[~, a, b, c, d, e, p0, p1, q0, q1, r0, r1] = getCoefficients(scheme, nu, D, h, tau);

%% Setup matrices for the scheme
[U_next, U_now, U_prev, F_next, F_now, F_prev] = getMatrices(Nx, a, b, c, d, e, p0, p1, q0, q1, r0, r1);

%% Integration loop
for k = 3:Nt
    % This part will implement the actual time-stepping scheme.
    rhs = -U_now*u(k-1, :).' - U_prev*u(k-2, :).' + ...
        F_next*f(k, :).' + F_now*f(k-1, :).' + F_prev*f(k-2, :).';
    u(k, :) = U_next \ rhs;
end

%% Apply boundary conditions if needed
u = [u, u(:, 1)];
end
