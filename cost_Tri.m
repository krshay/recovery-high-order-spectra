function [f] = cost_Tri(z, y, A, k1k2k3k4_map)
%cost_Tri Calculates the cost and grad for the optimization
%using the Trispectrum

N = size(z, 1) / 2;

T = calcTrispectrum(z(1:N) + 1j * z(N+1:end), k1k2k3k4_map);

f = (1 / 2) * sum(abs(y - A * T(:)).^2);
end