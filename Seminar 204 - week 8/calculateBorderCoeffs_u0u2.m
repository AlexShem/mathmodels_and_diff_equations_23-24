%% Declare symbolic variables
syms x t u(t, x) f(t, x)
syms D C h tau positive
syms alpha beta [2 4]
syms nu mu

alpha_mat = alpha;
alpha_mat(2, [1:1, 4:end]) = 0;
% alpha_mat(3, [1:1, 2:end]) = 0;
alpha = nonzeros(alpha_mat);

beta_mat = beta;
beta_mat(1, 3:end) = 0;
beta_mat(2, :) = 0;
beta = nonzeros(beta_mat);

%% Define test functions
test_funs = [x, x^3, x^4,...
    t*x, t*x^3,...
    t^2*x];

% h-tau offset grid
[h_mat, tau_mat] = meshgrid( ...
    (0 : size(alpha_mat, 2)-1)*h, ...
    (size(alpha_mat, 1)-1 : -1 : 0)*tau ...
    );

% Rod differential operator
P = @(U) diff(U, t, 2) - D*diff(diff(U, x, 2), t, 2) + C*diff(U, x, 4);

%% Build 2 systems
coeffs = sym('Eqs', [length(test_funs) + 2, 1]);

for k = 1 : numel(test_funs)
    u(t, x) = test_funs(k);
    f(t, x) = P(u);

    u_compact = sum(alpha_mat .* u(tau_mat, h_mat), 'all');
    f_compact = sum(beta_mat .* f(tau_mat, h_mat), 'all');

    coeffs(k, 1) = u_compact == f_compact;
end

%% Compute the coeffs of compact border
% Create two systems
systemBorder = coeffs;
systemPreBorder = coeffs;

% Introduce the normalisation
% Border equals to 1
systemBorder(end-1) = alpha1_1 == 1;
systemBorder(end) = alpha1_2 == 0;
% Pre-border equals to 1
systemPreBorder(end-1) = alpha1_1 == 0;
systemPreBorder(end) = alpha1_2 == 1;

% Solutions
coeffsBorder = solve(systemBorder, [alpha; beta]);
coeffsPreBorder = solve(systemPreBorder, [alpha; beta]);

%% Save the results
filename = 'compactBorder_u0u2.mat';
save(filename, "coeffsBorder", "coeffsPreBorder");
