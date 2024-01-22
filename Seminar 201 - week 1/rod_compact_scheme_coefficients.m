%% 5-5-5 Compact Scheme for the Rod Equation
clear
clc
syms a b c d e % Coeffs of the left side (approximating u)
syms p q r [1 2] % Coeffs of the right side (approximating f)
syms C D h tau nu mu positive % Positive parameters
syms x t u(t, x) f(t, x)

coefs_u_app = [a b c d e];
coefs_f_app = [p q r];

x_test = x.^[0 2 4 6];
t_test = t.^[0 2 4];
test_funs = transpose(x_test.' * t_test);
test_funs = test_funs(:);
% Keep 11 test functions
test_funs = test_funs(1:end-1);

%% Differential equation: Rod
du_fun = @(U) diff(U, t, 2) ...
    - D * diff(diff(U, x, 2), t, 2) ...
    + C * diff(U, x, 4);

coef_eqs = sym('Eqs', [length(test_funs), 1]);

for k = 1 : length(test_funs)
    u(t, x) = test_funs(k);
    f(t, x) = du_fun(u(t, x));
    
    u_comact = 1*(u(tau, 0) + u(-tau, 0)) + ...
        b * u(0, 0) + ...
        c * (u(0, h) + u(0, -h)) + ...
        a * (u(tau, h) + u(tau, -h) + u(-tau, h) + u(-tau, -h)) + ...
        d * (u(0, 2*h) + u(0, -2*h)) + ...
        e * (u(tau, 2*h) + u(tau, -2*h) + u(-tau, 2*h) + u(-tau, -2*h));
    f_comact = r1 * f(0, 0) + ...
        r2 * (f(tau, 0) + f(-tau, 0)) + ...
        q1 * (f(0, h) + f(0, -h)) + ...
        q2 * (f(tau, h) + f(tau, -h) + f(-tau, h) + f(-tau, -h)) + ...
        p1 * (f(0, 2*h) + f(0, -2*h)) + ...
        p2 * (f(tau, 2*h) + f(tau, -2*h) + f(-tau, 2*h) + f(-tau, -2*h));
    
    coef_eqs(k, 1) = u_comact == f_comact;
end

%% Solution
comp_eqs = solve(coef_eqs, [coefs_u_app, coefs_f_app]);

cfs = struct2array(comp_eqs).';
cfs = simplify(subs(cfs, [C*tau^2, D], [nu*h^4, mu*h^2]));

a_st = 72*(41 + 30*nu + 90*mu);
cfs_norm = simplify(cfs * a_st);
