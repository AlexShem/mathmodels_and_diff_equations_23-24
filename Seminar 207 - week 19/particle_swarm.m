% Define the fitness function
fitnessFunction = @(x) x(1)^2 + x(2)^2; % Example: minimize the sphere function
fitnessFunction = @(x) exp(-(x(1).^2+x(2).^2)/10).*sin(x(1)).*cos(x(2));

% Parameters
numParticles = 30;  % Number of particles
numDimensions = 2;  % Number of dimensions
maxIterations = 100; % Maximum number of iterations
w = 0.5;            % Inertia weight
c1 = 1.5;           % Cognitive (personal) coefficient
c2 = 1.5;           % Social (global) coefficient

% Initialize particle positions and velocities
positions = rand(numParticles, numDimensions) * 10 - 5; % Random initial positions in range [-5, 5]
velocities = rand(numParticles, numDimensions) * 2 - 1; % Random initial velocities in range [-1, 1]
personalBestPositions = positions; % Initial personal best positions
personalBestScores = arrayfun(@(idx) fitnessFunction(positions(idx, :)), 1:numParticles);
[globalBestScore, bestIndex] = min(personalBestScores);
globalBestPosition = personalBestPositions(bestIndex, :);

% PSO main loop
for iter = 1:maxIterations
    % Update velocities and positions
    r1 = rand(numParticles, numDimensions);
    r2 = rand(numParticles, numDimensions);
    velocities = w * velocities + c1 * r1 .* (personalBestPositions - positions) + c2 * r2 .* (globalBestPosition - positions);
    positions = positions + velocities;
    
    % Evaluate fitness
    scores = arrayfun(@(idx) fitnessFunction(positions(idx, :)), 1:numParticles);
    
    % Update personal best positions and scores
    for i = 1:numParticles
        if scores(i) < personalBestScores(i)
            personalBestScores(i) = scores(i);
            personalBestPositions(i, :) = positions(i, :);
        end
    end
    
    % Update global best position and score
    [newBestScore, bestIndex] = min(personalBestScores);
    if newBestScore < globalBestScore
        globalBestScore = newBestScore;
        globalBestPosition = personalBestPositions(bestIndex, :);
    end
    
    % Display the iteration number and the current best score
    disp(['Iteration ' num2str(iter) ': Best Score = ' num2str(globalBestScore)]);
end

% Display the results
disp('Optimization completed.');
disp(['Best position: ' num2str(globalBestPosition)]);
disp(['Best score: ' num2str(globalBestScore)]);
