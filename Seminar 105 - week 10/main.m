%% Main parameters
D = 0.1; % Diffusion coefficient
L = 1; % Length of the rod
T = 1.5; % Time of integration
Nx = 50; % Number of spatial point
nu = 0.5; % Curant parameter

% nu = D * tau / h^2

%% Secondary parameters
x = linspace(0, L, Nx); % Spatial mesh
h = x(2) - x(1); % Spatial step

tau = nu * h^2 / D; % Temportal step
Nt = ceil(T/tau) + 1; % Number of temporal steps
t = 0:tau:((Nt-1)*tau);

%% Initail condition
u_0 = exp(-(x - .5).^2 ./ .05);
% u_0 = sin(2*pi*x/L);

%% Transition matrices
U_now = (1 - 2*nu) * eye(Nx) + ...
    diag(nu*ones(Nx-1, 1), 1) + ...
    diag(nu*ones(Nx-1, 1), -1);
U_next = eye(Nx);

%% Boundary conditions
% Dirichlet conditions u(t, 0) = u(t, L) = 0
U_now(1, :) = 0;
U_now(end, :) = 0;

% % Neuman conditions: d_x u (t, 0) = d_x u (t, L) = 0
% % u(t, 0) = u(t, h)
% % u(t, L) = u(t, L-h)
% U_now(1, :) = 0;
% U_now(end, :) = 0;
% 
% U_next(1, 2) = -1;
% U_next(end, end-1) = -1;

%% Inegration of the system
u = zeros(Nt, Nx);
u(1, :) = u_0;

for k = 2 : Nt
    u(k, :) = U_next \ U_now * u(k-1, :).';
end

%% Visualisation
[x_mesh, t_mesh] = meshgrid(x, t);
figure(1)
surf(x_mesh, t_mesh, u, 'LineStyle','none')
xlabel('x');
ylabel('t');

%% Animation (lazy)
figure(2)
for k = [1 : floor(Nt/100) : Nt, Nt]
    plot(x, u(k, :));
    xlabel('x');
    title(['t = ', num2str(tau*(k-1))]);
    axis([0, L, -1, 1]);
    drawnow
end
