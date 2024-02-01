function isNewFigure = plotErrorNorms(Nx_values, C_norms, L2_norms, figureId)
% Initialize the flag
isNewFigure = true;

% Check if figureId is provided
if nargin >= 4 && ~isempty(findall(0, 'Type', 'Figure', 'Number', figureId))
    figure(figureId); % Make the figure with figureId the current figure
    isNewFigure = false; % The figure already exists
else
    % Open a new figure window and use the provided figureId if it exists
    if nargin >= 4
        figure(figureId);
    else
        figure; % Open a new figure without a specific id
    end
end

% Create subplots
% Subplot for C norms
subplot(1, 2, 1);
if ~isNewFigure
    hold on;
end
loglog(Nx_values, C_norms, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Number of spatial points, $N_x$', 'Interpreter', 'latex');
ylabel('C Norm', 'Interpreter', 'latex');
title('C Norm vs. Number of Spatial Points', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);
hold off;

% Subplot for L2 norms
subplot(1, 2, 2);
if ~isNewFigure
    hold on;
end
loglog(Nx_values, L2_norms, 's-', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Number of spatial points, $N_x$', 'Interpreter', 'latex');
ylabel('$L^2$ Norm', 'Interpreter', 'latex');
title('$L^2$ Norm vs. Number of Spatial Points', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);
hold off;

% Enhance the figure layout with a super title
sgtitle('Error Norms vs. Number of Spatial Points', 'Interpreter', 'latex');

% Reset hold state to off for future plots
% hold off;

% Optionally return isNewFigure flag
if nargout == 0
    clear isNewFigure; % If not capturing the output, don't display it
end
end
