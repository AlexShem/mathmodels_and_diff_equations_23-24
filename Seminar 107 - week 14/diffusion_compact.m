%% Main parameters
D = .1;
L = 2*pi;
T = 10;

Nx = 100;
nu = .25; % nu = D * tau / h^2

%% Secondary parameters
x = linspace(0, L, Nx + 1);
h = x(2) - x(1);

tau = nu * h^2 / D;
Nt = ceil(T/tau) + 1;
t = linspace(0, T, Nt).';

u_true = @(t, x) sin(t) * cos(x);
u_0 = u_true(0, x);
f = @(t, x) (D*sin(t) + cos(t))*cos(x);

%% Transition matrices
alpha = 4*(6*nu + 5);
a = 2*(1 - 6*nu)/alpha;
b = -4*(5 - 6*nu)/alpha;
c = -2*(1 + 6*nu)/alpha;
p1 = tau/alpha;
p2 = tau/alpha;
q1 = 10*tau/alpha;
q2 = 10*tau/alpha;

U_now = b*eye(Nx) + ...
    diag(c * ones(Nx-1, 1), 1) + ...
    diag(c * ones(Nx-1, 1),-1);
U_next = eye(Nx) + ...
    diag(a * ones(Nx-1, 1), 1) + ...
    diag(a * ones(Nx-1, 1),-1);
f_next = diag(q2*ones(Nx, 1)) + ...
    diag(p2*ones(Nx-1, 1), 1) + diag(p2*ones(Nx-1, 1), -1);
f_now = diag(q1*ones(Nx, 1)) + ...
    diag(p1*ones(Nx-1, 1), 1) + diag(p1*ones(Nx-1, 1), -1);

%% Periodic border condition
U_next(1, end) = a;
U_next(end, 1) = a;
U_now(1, end) = c;
U_now(end, 1) = c;
f_next(1, end) = p2;
f_next(end, 1) = p2;
f_now(1, end) = p1;
f_now(end, 1) = p1;

%% Integration
U = zeros(Nt, Nx);
U(1, :) = u_0(1:end-1);
f_val = f(t, x(1:end-1));

for k = 2 : Nt
    right_part = -U_now*U(k-1, :).' + ...
        f_next*f_val(k, :).' + f_now*f_val(k-1, :).';
    U(k, :) = U_next \ right_part;
end
U = [U, U(:, 1)];

%% Visualisation
figure(2)

for k = [1 : floor(Nt/200) : Nt, Nt]
    plot(x, U(k, :));
    hold on;
    fplot(@(x)u_true(t(k), x), [0, L]);
    hold off;
    axis([0 L -1 1]);
    title(['t = ', num2str((k-1)*tau)]);
    drawnow;
end
