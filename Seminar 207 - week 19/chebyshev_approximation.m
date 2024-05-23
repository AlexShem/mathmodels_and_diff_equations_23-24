% Step 1: Define the function handle
f = @(x) exp(x); % Example function
f = @(x) x.^3; % Example function

% Step 2: Generate Chebyshev nodes
poly_degree = 1; % Degree of the polynomial
n = poly_degree + 1;
k = 0:n-1;
x_nodes = cos((2*k+1)*pi/(2*n));

% Step 3: Evaluate the function at the Chebyshev nodes
y_values = f(x_nodes);

% Step 4: Compute Chebyshev coefficients
a = zeros(1, n);
for j = 1:n
    a(j) = (2/n) * sum(y_values .* cos((j-1) * acos(x_nodes)));
end
a(1) = a(1) / 2; % Adjust the constant term

% Step 5: Construct the interpolating polynomial
% Define the Chebyshev polynomials
T = @(n, x) cos(n * acos(x)); 

% Evaluate the interpolating polynomial at points in [-1, 1]
x = linspace(-1, 1, 1000);
P_n = a(1) * ones(size(x)); % Start with the constant term
for j = 2:n
    P_n = P_n + a(j) * T(j-1, x);
end

% Plot the original function and the interpolating polynomial
figure(1);
plot(x, f(x), 'b', 'DisplayName', 'Original Function');
hold on;
plot(x, P_n, 'r--', 'DisplayName', 'Chebyshev Polynomial Approximation');
legend show;
title('Chebyshev Polynomial Approximation');
xlabel('x');
ylabel('f(x) and P_n(x)');
grid on;
hold off;
