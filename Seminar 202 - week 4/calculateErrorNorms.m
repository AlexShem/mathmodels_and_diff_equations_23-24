function [C_norm, L2_norm] = calculateErrorNorms(u, u_ref, x, t, T)
% u - matrix of the solution at discrete times t and space x
% u_ref - function handle @(t, x) of the reference solution
% x - spatial grid
% t - array of times at which the solution u is calculated
% T - target times for error calculation; can be 'last' or a non-negative numeric vector

% Initialize variables to store maximum (C) and L2 norms for each requested time T
if ischar(T) && strcmp(T, 'last')
    T = t(end); % Use the last time step
end

% Ensure T is a row vector for consistent processing
if iscolumn(T)
    T = T';
end

% Preallocate arrays for norms
C_norm = zeros(size(T));
L2_norm = zeros(size(T));

% Spatial step for L2 norm calculation
dx = x(2) - x(1);

% Loop over each requested time in T
for i = 1:length(T)
    targetTime = T(i);

    % Find the closest times in t before and after targetTime
    if targetTime >= t(end)
        % If targetTime is beyond the last calculated time, use the last time step
        u_numerical = u(end, :);
    else
        % Interpolate u at targetTime
        u_numerical = interp2(x, t, u, x, targetTime, "spline");
    end

    % Calculate the reference solution at targetTime
    u_reference = u_ref(targetTime, x);

    % Calculate error
    error = abs(u_numerical - u_reference);

    % Calculate norms
    C_norm(i) = max(error);
    L2_norm(i) = sqrt(sum((error.^2) * dx));
end
end
