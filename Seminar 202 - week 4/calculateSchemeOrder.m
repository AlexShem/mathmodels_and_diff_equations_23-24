%% Main Script to Calculate Scheme Order
% This script calculates the order of a computational scheme (CN or "555")
% by plotting log-log plots of C and L^2 norms of the difference between
% the obtained solution u and the theoretical solution u_ref depending on
% the number of spatial steps.

clear; clc;
% close all;

%% Primary Parameters
rho = 7900; % Density of steel (kg/m^3)
R = 10^-2; % Radius of the rod (m)
E = 210e9; % Young's modulus (Pa)
L = 2*pi; % Length of the rod (m)
T = 1; % Total simulation time (s)
nu = 0.1; % Constant parameter

scheme = "CN"; % Computational scheme ("CN" or "555")
% scheme = "555"; % Computational scheme ("CN" or "555")
problemId = 1; % Identificator of the reference solution (1 or 2)

%% Spatial discretizations to analyze
Nx_values = unique(round(logspace(1, 2, 10)));

%% Initialize arrays for error norms
C_norms = zeros(size(Nx_values));
L2_norms = zeros(size(Nx_values));

%% Loop over spatial discretizations
for i = 1:length(Nx_values)
    Nx = Nx_values(i);
    [C_norm, L2_norm] = performAnalysis(Nx, L, T, rho, R, E, nu, scheme, problemId);
    C_norms(i) = C_norm;
    L2_norms(i) = L2_norm;
end

%% Visualization of Error Norms
plotErrorNorms(Nx_values, C_norms, L2_norms, 2);
