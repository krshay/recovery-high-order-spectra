function [f] = cost_Bispectrum(z, y, A, k1k2k3_map)
%cost_Bispectrum Calculates the cost function for the optimization
%using the Bispectrum; eq. (4.2)

N = size(z, 1);

B = calcBispectrum(z, k1k2k3_map);

B_flat = reshape(B, N^2, 1);

f = (1 / 2) * sum(abs(y - A * B_flat).^2);
end