function [f] = cost_Bi(z, y, A, k1k2k3_map)
%cost_Bispectrum Calculates the cost function for the optimization
%using the Bispectrum; eq. (4.2)

N = size(z, 1) / 2;

B = calcBispectrum(z(1:N) + 1j * z(N+1:end), k1k2k3_map);

f = (1 / 2) * sum(abs(y - A * B(:)).^2);
end