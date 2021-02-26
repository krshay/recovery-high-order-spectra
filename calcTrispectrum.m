function T = calcTrispectrum(x, k1k2k3k4_map)
%calcTrispectrum Calculates the Trispectrum; eq. (1.3)

%   x is a complex N*1 vector of Fourier coefficients
%   k1k2k3k4_map is a N*N matrix containing the quadruple of frequencies

%   B = x[k1]x[k2]x[k3]x[-k1-k2-k3]

N = size(x, 1);

T = zeros(N, N, N);

for i=1:size(k1k2k3k4_map, 1)
    k1 = k1k2k3k4_map(i, 1);
    k2 = k1k2k3k4_map(i, 2);
    k3 = k1k2k3k4_map(i, 3);
    k4 = k1k2k3k4_map(i, 4);
    T(k1+1, k2+1, k3+1) = x(k1+1) * x(k2+1) * x(k3+1) * x(k4+1);
end
end