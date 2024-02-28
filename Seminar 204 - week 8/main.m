%% Primary parameters
% L - length of the rod;
% T - desired final time;
% rho, R, E - parameters of the rod;
% nu - dimensionless parameter, nu = (E*R^2/rho) * tau^2 / h^4;
% scheme - one of the computational schemes, e.g., "CN" or "555";
% omega - time frequency
% Nx - number of spatial points;
L = pi;
T = 1;
rho = 7900;
R = 1e-2;
E = 210e9;
nu = 0.1;
omega = pi;
Nx = 45;

%% Scheme and Boundary Condition Choice
% scheme = "CN";
scheme = "555";
% boundary_conditions = "u0u1";
boundary_conditions = "u0u2";
% approximation_type = "naive";
approximation_type = "compact";

%% Secondary parameters
% u_ref, f_ref - function handles to the reference solution and the external force.
D = R^2;
C = E*R^2/rho;
x = linspace(0, L, Nx);
h = x(2) - x(1);
mu = D / h^2;
tau = sqrt(nu*h^4/C);
Nt = ceil(T/tau) + 1;
t = (0 : tau : tau*(Nt-1)).';

if boundary_conditions == "u0u1"
    c_t = @(t) cos(omega * t);
    u_ref = @(t, x) c_t(t) .* sin(x).^2;
    f_ref = @(t, x) cos(omega.*t).*((omega.^2.*cos(x.*2.0))./2.0-omega.^2./2.0-C.*cos(x.*2.0).*8.0+D.*omega.^2.*cos(x.*2.0).*2.0);
elseif boundary_conditions == "u0u2"
    c_t = @(t) cos(omega * t);
    u_ref = @(t, x) c_t(t) .* sin(x);
    f_ref = @(t, x) -cos(omega.*t).*sin(x).*(-C+D.*omega.^2+omega.^2);
end

%% Coefficients of the computational scheme
[~, a, b, c, d, e, p0, p1, q0, q1, r0, r1] = getCoefficients(scheme, nu, D, h, tau);

%% Construct computational matrices
% Here, getMatrixes() applies boundary conditions
[U_next, U_now, U_prev, F_next, F_now, F_prev] = getMatrices(Nx, nu, mu, C, D, h, tau, boundary_conditions, approximation_type, a, b, c, d, e, p0, p1, q0, q1, r0, r1);

%% Initial conditions and External force
u = zeros(Nt, Nx); % Initialize solution matrix
u(1, :) = u_ref(0, x); % Apply initial condition at t=0
u(2, :) = u_ref(tau, x); % Apply initial condition at t=tau

[x_gr, t_gr] = meshgrid(x, t);
u_star = u_ref(t_gr, x_gr);
f = f_ref(t_gr, x_gr);

%% Integration loop
for k = 3:Nt
    % This part will implement the actual time-stepping scheme.
    rhs = -U_now*u(k-1, :).' - U_prev*u(k-2, :).' + ...
        F_next*f(k, :).' + F_now*f(k-1, :).' + F_prev*f(k-2, :).';
    u(k, :) = U_next \ rhs;
end

%% Animation
figure(1)
for k = [1 : 100 : Nt, Nt]
    plot(x, u(k, :), '-b');
    hold on;
    plot(x, u_star(k, :), '--r');
    hold off;
    xlabel('$x$', FontSize=16, Interpreter="latex");
    title(['$t = $', num2str(round(t(k), 3))], FontSize=16, Interpreter="latex");
    axis([0, L, -1.1, 1.1]);
    drawnow limitrate nocallbacks
end

%% Error
C_norm_t = max(abs(u - u_star), [], 2);
figure(3);
semilogy(t, C_norm_t);
hold on;
