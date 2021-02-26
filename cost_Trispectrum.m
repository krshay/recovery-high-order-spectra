function [f] = cost_Trispectrum(z, y, A, k1k2k3k4_map)
%func_cost_Trispectrum Calculates the cost and grad for the optimization
%using the Trispectrum

N = size(z, 1);

T = calcTrispectrum(z, k1k2k3k4_map);

T_flat = reshape(T, N^3, 1);

f = (1 / 2) * sum(abs(y - A * T_flat).^2);
end