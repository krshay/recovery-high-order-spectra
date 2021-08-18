function [z, cst] = optimize(z_init, y, A, map)
% Requires Manopt: http://www.manopt.org/

if size(map, 2) == 3
    fun = @(z) cost_Bi(z, y, A, map);
else
    fun = @(z) cost_Tri(z, y, A, map);
end
options = optimoptions(@fminunc,'Display','off','OptimalityTolerance', 1e-15, 'StepTolerance', 1e-15, 'MaxFunctionEvaluations', 10000);

[z, cst] =  fminunc(fun, z_init, options);
end