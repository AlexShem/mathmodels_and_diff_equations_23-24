function [alpha, a, b, c, d, e, p0, p1, q0, q1, r0, r1] = getCoefficients(scheme, nu, D, h, tau)
%% Initialize coefficients
alpha = 0; a = 0; b = 0; c = 0; d = 0; e = 0;
p0 = 0; p1 = 0; q0 = 0; q1 = 0; r0 = 0; r1 = 0;

%% Calculate mu and other scheme-specific parameters
mu = D / h^2;

%% Coefficients based on the scheme
if strcmp(scheme, "CN")
    % Crank-Nicolson Scheme Coefficients
    alpha = 1 + 3*nu + 2*mu;
    a = -(2*nu + mu) / alpha;
    b = -(2 + 4*mu) / alpha;
    c = 2*mu / alpha;
    d = 0; % Not used in CN scheme, set to 0 for clarity
    e = nu/2 / alpha;
    % Right part coefficients (not used in CN, set to 0 for consistency)
    p0 = 0; p1 = 0; q0 = 0; q1 = 0; r0 = 0; r1 = tau^2 / 2 / alpha;
elseif strcmp(scheme, "555")
    % "5-5-5" Compact Scheme Coefficients
    alpha = 72*(41 + 30*nu + 90*mu);
    a = 96*(7 - 15*nu - 30*mu) / alpha;
    b = -144*(41 - 150*nu + 90*mu) / alpha;
    c = -192*(7 + 75*nu - 30*mu) / alpha;
    d = -24*(1 - 150*nu - 30*mu) / alpha;
    e = 12*(1 + 30*nu - 30*mu) / alpha;
    % Right part coefficients
    p0 = 10*tau^2 / alpha;
    q0 = 560*tau^2 / alpha;
    r0 = 2460*tau^2 / alpha;
    p1 = tau^2 / alpha;
    q1 = 56*tau^2 / alpha;
    r1 = 246*tau^2 / alpha;
else
    error('Unsupported scheme: %s', scheme);
end
end
