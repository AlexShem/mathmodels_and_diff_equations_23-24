% GRAD_DESCENT Performs gradient descent to minimize a function.
% [xmin, fmin, niter, path] = grad_descent(x0, f, tol, maxiter) performs
% gradient descent on the function f, starting from the initial guess x0,
% until the change in f(x) is less than the tolerance tol or the maximum
% number of iterations maxiter is reached.
%
% Inputs:
% x0 - Initial guess for the minimum.
% f - Function to be minimized.
% tol - Tolerance for convergence (default: 1e-8).
% maxiter - Maximum number of iterations (default: 50).
%
% Outputs:
% xmin - Minimum x.
% fmin - Minimum value of f.
% niter - Number of iterations.
% path - Path of x at each iteration.

function [xmin, fmin, niter, path] = grad_descent(x0, f, tol, maxiter)
if nargin < 4
    maxiter = 50; % Default maximum number of iterations
end
if nargin < 3 || (nargin == 4 && isempty(tol))
    tol = 1e-8; % Default tolerance
end

x0 = x0(:); % Ensure x0 is a column vector
path = zeros(length(x0), maxiter + 1); % Initialize path
path(:, 1) = x0; % Store initial guess in path

lam = 1; % Initial step size
grad_x = grad_f(x0, f); % Calculate gradient at x0
x = x0 - grad_x; % Update x
iter = 1; % Initialize iteration counter
path(:, iter + 1) = x; % Store x in path

df = abs(f(x) - f(x0)); % Calculate change in f(x)

% Main loop
while (df > tol) && iter < maxiter
    grad_x = grad_f(x, f); % Calculate gradient at x

    x0 = x; % Update x0
    x = x0 - lam * grad_x; % Update x
    % If f(x) is not decreasing, halve the step size
    while f(x) > f(x0)
        lam = lam/2;
        x = x0 - lam * grad_x;
    end
    iter = iter + 1; % Increment iteration counter
    path(:, iter + 1) = x; % Store x in path

    df = abs(f(x) - f(x0)); % Calculate change in f(x)
end

% If maximum number of iterations reached without convergence, issue warning
if iter == maxiter
    warning('Maximum number of iterations reached without convergence.')
end

xmin = x; % Minimum x
fmin = f(x); % Minimum value of f
niter = iter; % Number of iterations
path = path(:, 1 : niter + 1);

function grad = grad_f(x, f)
d = length(x);
h = 1e-8;
grad = arrayfun(@(ind) (f(x + h*((1:d).' == ind)) - f(x - h*((1:d).' == ind)))/(2*h), (1:d).');
grad = grad / norm(grad, 2);
