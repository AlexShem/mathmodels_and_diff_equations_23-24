function [U_next, U_now, U_prev, F_next, F_now, F_prev] = getMatrices(Nx, a, b, c, d, e, p0, p1, q0, q1, r0, r1)
U_next = eye(Nx) + diag(a*ones(Nx-1, 1), 1) + diag(a*ones(Nx-1, 1), -1) + ...
    diag(e*ones(Nx-2, 1), 2) + diag(e*ones(Nx-2, 1), -2);
U_now = b*eye(Nx) + diag(c*ones(Nx-1, 1), 1) + diag(c*ones(Nx-1, 1), -1) + ...
    diag(d*ones(Nx-2, 1), 2) + diag(d*ones(Nx-2, 1), -2);

F_next = r1*eye(Nx) + diag(q1*ones(Nx-1, 1), 1) + diag(q1*ones(Nx-1, 1), -1) + ...
    diag(p1*ones(Nx-2, 1), 2) + diag(p1*ones(Nx-2, 1), -2);
F_now = r0*eye(Nx) + diag(q0*ones(Nx-1, 1), 1) + diag(q0*ones(Nx-1, 1), -1) + ...
    diag(p0*ones(Nx-2, 1), 2) + diag(p0*ones(Nx-2, 1), -2);

% Periodic boundary conditions
U_next(1, end-1:end) = [e a];
U_next(2, end) = e;
U_next(end-1, 1) = e;
U_next(end, 1:2) = [a e];

U_now(1, end-1:end) = [d c];
U_now(2, end) = d;
U_now(end-1, 1) = d;
U_now(end, 1:2) = [c d];

F_next(1, end-1:end) = [p1 q1];
F_next(2, end) = p1;
F_next(end-1, 1) = p1;
F_next(end, 1:2) = [q1 p1];

F_now(1, end-1:end) = [p0 q0];
F_now(2, end) = p0;
F_now(end-1, 1) = p0;
F_now(end, 1:2) = [q0 p0];

U_prev = U_next;
F_prev = F_next;
end
