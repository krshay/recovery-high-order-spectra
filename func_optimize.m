function [z, problem, cst] = func_optimize(z_init, y, A, map)
% Requires Manopt: http://www.manopt.org/
N = size(z_init, 1);

problem.M = euclideancomplexfactory(N, 1);

if size(map, 2) == 3
    problem.cost = @(z) cost_Bispectrum(z, y, A, map);
else
    problem.cost = @(z) cost_Trispectrum(z, y, A, map);
end
options.tolgradnorm = 1e-4;
options.verbosity = 0;
options.maxiter = 100000;
options.minstepsize = 1e-10;
[z, cst] = steepestdescent(problem, z_init, options);
end