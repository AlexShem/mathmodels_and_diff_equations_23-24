function [u_ref, f_ref] = getReferenceSolution(problemId, C, D)
if problemId == 1
    u_ref = @(t, x) sin(x + t*sqrt(C/(D+1))); % Reference solution
    f_ref = @(t, x) zeros(size(x));
else
    error('Unsupported problem ID: %d', problemId);
end
end
