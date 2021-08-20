function B = calcBispectrum(x, k1k2k3_map)
%calcBispectrum Calculates the Bispectrum; eq. (1.1)

%   x is a complex N*1 vector of Fourier coefficients
%   k1k2k3_map is a N^2*3 matrix containing the triplets of frequencies

%   B = x[k1]x[k2]x[-k1-k2]

N = size(x, 1);

B = reshape(x(k1k2k3_map(:, 1) + 1) .* x(k1k2k3_map(:, 2) + 1) .* x(k1k2k3_map(:, 3) + 1), N, N);
end
