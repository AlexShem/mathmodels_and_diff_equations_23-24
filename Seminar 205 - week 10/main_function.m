function main_function(L, T, tau, Nx, comp_correction, scheme)
    % Set up parameters
    [x, h, Nt, t, u_0, f, ax] = setup_parameters(L, T, tau, Nx);

    % Perform integration
    U = perform_integration(Nt, Nx, u_0, f, tau, h, comp_correction, scheme);

    % Visualize results
    visualize_results(U, Nt, tau, x, ax);
end

function [x, h, Nt, t, u_0, f, ax] = setup_parameters(L, T, tau, Nx)
    x = linspace(0, L, Nx + 1);
    h = x(2) - x(1);
    Nt = ceil(T/tau) + 1;
    T = tau*(Nt - 1);
    t = 0 : tau : T;
    u_0 = sin(x);
    ax = [0 L -2 2];
    f = @(u) u.^2/2;
end

function U = perform_integration(Nt, Nx, u_0, f, tau, h, comp_correction, scheme)
    U = zeros(Nt, Nx);
    U(1, :) = u_0(1 : end-1);
    for k = 2 : Nt
        u = U(k-1, :);
        switch scheme
            case 'explicit_euler'
                u_ex = explicit_euler(u, f, tau, h);
            case 'maccormack'
                u_ex = maccormack(u, f, tau, h);
            case 'lax_wendroff'
                u_ex = lax_wendroff(u, f, tau, h);
            otherwise
                error('Invalid scheme');
        end
        if comp_correction
            comp_eps = compact_correction(u, u_ex, h, tau);
        else
            comp_eps = zeros(length(u), 1);
        end
        U(k, :) = u_ex + comp_eps.';
    end
    U = [U, U(:, 1)];
end

function visualize_results(U, Nt, tau, x, ax)
    figure(2)
    for k = [1 : floor(Nt/200) : Nt, Nt]
        plot(x, U(k, :));
        axis(ax);
        title(['t = ', num2str((k-1)*tau)]);
        drawnow;
    end
end
