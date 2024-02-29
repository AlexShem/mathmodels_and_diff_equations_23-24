function [U_next, U_now, U_prev, F_next, F_now, F_prev] = getMatrices(Nx, nu, mu, C, D, h, tau, boundary_conditions, approximation_type, a, b, c, d, e, p0, p1, q0, q1, r0, r1)
U_next = eye(Nx) + diag(a*ones(Nx-1, 1), 1) + diag(a*ones(Nx-1, 1), -1) + ...
    diag(e*ones(Nx-2, 1), 2) + diag(e*ones(Nx-2, 1), -2);
U_now = b*eye(Nx) + diag(c*ones(Nx-1, 1), 1) + diag(c*ones(Nx-1, 1), -1) + ...
    diag(d*ones(Nx-2, 1), 2) + diag(d*ones(Nx-2, 1), -2);
U_prev = eye(Nx) + diag(a*ones(Nx-1, 1), 1) + diag(a*ones(Nx-1, 1), -1) + ...
    diag(e*ones(Nx-2, 1), 2) + diag(e*ones(Nx-2, 1), -2);

F_next = r1*eye(Nx) + diag(q1*ones(Nx-1, 1), 1) + diag(q1*ones(Nx-1, 1), -1) + ...
    diag(p1*ones(Nx-2, 1), 2) + diag(p1*ones(Nx-2, 1), -2);
F_now = r0*eye(Nx) + diag(q0*ones(Nx-1, 1), 1) + diag(q0*ones(Nx-1, 1), -1) + ...
    diag(p0*ones(Nx-2, 1), 2) + diag(p0*ones(Nx-2, 1), -2);
F_prev = r1*eye(Nx) + diag(q1*ones(Nx-1, 1), 1) + diag(q1*ones(Nx-1, 1), -1) + ...
    diag(p1*ones(Nx-2, 1), 2) + diag(p1*ones(Nx-2, 1), -2);

if approximation_type == "compact"
    if boundary_conditions == "u0u1"
        alpha_mat_border = [1, 0, 0, 0;
            0, 0, 0, 0];
        alpha_mat_preborder = [0, 1, -1/2, 1/9;
            0, 0, 0, 0];
        beta_mat_border = [0, 0, 0, 0;
            0, 0, 0, 0];
        beta_mat_preborder = [-h^4*(2*mu-1)/(12*C), -h^4*mu/(6*C), 0, 0;
            0, 0, 0, 0];
    elseif boundary_conditions == "u0u2"
        alpha_mat_border = [1, 0, 0, 0;
            0, 0, 0, 0];
        alpha_mat_preborder = [0, 1, -4/5, 1/5;
            0, 0, 0, 0];
        beta_mat_border = [0, 0, 0, 0;
            0, 0, 0, 0];
        beta_mat_preborder = [11*h^4/(60*C), 0, 0, 0;
            0, 0, 0, 0];
    end
elseif approximation_type == "naive"
    if boundary_conditions == "u0u1"
        alpha_mat_border = [1, 0, 0, 0;
            0, 0, 0, 0];
        alpha_mat_preborder = [-1, 1, 0, 0;
            0, 0, 0, 0];
        beta_mat_border = [0, 0, 0, 0;
            0, 0, 0, 0];
        beta_mat_preborder = [0, 0, 0, 0;
            0, 0, 0, 0];
    elseif boundary_conditions == "u0u2"
        alpha_mat_border = [1, 0, 0, 0;
            0, 0, 0, 0];
        alpha_mat_preborder = [2, -5, 4, -1;
            0, 0, 0, 0];
        % alpha_mat_preborder = [1, -2, 1, 0;
        %     0, 0, 0, 0];
        beta_mat_border = [0, 0, 0, 0;
            0, 0, 0, 0];
        beta_mat_preborder = [0, 0, 0, 0;
            0, 0, 0, 0];
    end
end

U_next(1, 1:4) = alpha_mat_border(1, :);
U_next(2, 1:4) = alpha_mat_preborder(1, :);
U_next(end-1, end-3:end) = flip(alpha_mat_preborder(1, :));
U_next(end, end-3:end) = flip(alpha_mat_border(1, :));

U_now(1, 1:4) = alpha_mat_border(2, :);
U_now(2, 1:4) = alpha_mat_preborder(2, :);
U_now(end-1, end-3:end) = flip(alpha_mat_preborder(2, :));
U_now(end, end-3:end) = flip(alpha_mat_border(2, :));

U_prev([1:2, end-1:end], :) = 0;

F_next(1, 1:4) = beta_mat_border(1, :);
F_next(2, 1:4) = beta_mat_preborder(1, :);
F_next(end-1, end-3:end) = flip(beta_mat_preborder(1, :));
F_next(end, end-3:end) = flip(beta_mat_border(1, :));

F_now(1, 1:4) = beta_mat_border(2, :);
F_now(2, 1:4) = beta_mat_preborder(2, :);
F_now(end-1, end-3:end) = flip(beta_mat_preborder(2, :));
F_now(end, end-3:end) = flip(beta_mat_border(2, :));

F_prev([1:2, end-1:end], :) = 0;
end
