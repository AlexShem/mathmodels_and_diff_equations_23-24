function [C_norm, L2_norm] = performAnalysis(Nx, L, T, rho, R, E, nu, scheme, problemId)
% Define the reference solution function
D = R^2;
C = E*R^2/rho;
[u_ref, f_ref] = getReferenceSolution(problemId, C, D);

% Solve the rod equation for the given parameters and scheme
[t, u] = solveRodEquation(Nx, L, T, rho, R, E, nu, scheme, u_ref, f_ref);

% Generate the spatial grid
x = linspace(0, L, Nx + 1);

% Calculate error norms
% [C_norm, L2_norm] = calculateErrorNorms(u, u_ref, x, t, 'last');
[C_norm, L2_norm] = calculateErrorNorms(u, u_ref, x, t, T);
end
