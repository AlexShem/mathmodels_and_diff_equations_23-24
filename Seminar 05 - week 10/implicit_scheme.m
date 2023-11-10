%% Main parameters
D = .1;
L = 1;
T = 1.5;

Nx = 50;
nu = .25; % nu = D * tau / h^2

%% Secondary parameters
x = linspace(0, L, Nx);
h = x(2) - x(1);

tau = nu * h^2 / D;
Nt = ceil(T/tau) + 1;

u_0 = exp(-(x - .5).^2 ./ .05);

%% Transition matrices
U_now = -eye(Nx);
U_next = (1+2*nu)*eye(Nx) + ...
    diag(-nu * ones(Nx-1, 1), 1) + ...
    diag(-nu * ones(Nx-1, 1),-1);

%% Dirichlet border condition
U_next(1, :) = [1, zeros(1, Nx - 1)];
U_next(end, :) = [zeros(1, Nx - 1), 1];

U_now(1, :) = 0;
U_now(end, :) = 0;

%% Integration
U = zeros(Nt, Nx);
U(1, :) = u_0;

for k = 2 : Nt
    U(k, :) = -U_next \ U_now * U(k - 1, :).';
end

%% Visualisation
figure(2)

for k = [1 : floor(Nt/100) : Nt, Nt]
    plot(x, U(k, :));
    axis([0 L 0 1]);
    title(['t = ', num2str((k-1)*tau)]);
    drawnow;
end
