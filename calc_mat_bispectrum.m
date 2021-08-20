function A = calc_mat_bispectrum(K, N, k1k2k3_map)
%calc_mat_bispectrum Calculates the sensing matrix that consists
%of bispectra, for experiment 3a

%   K is the number of samples
%   N is the length of the signal
%   k1k2k3_map is a N^2*3 matrix containing the triplets of frequencies

A = zeros(K, N^2);
for k=1:K
    sig = randn(N, 1) + 1j * randn(N, 1);
    z = fft(sig);
    B_z = calcBispectrum(z, k1k2k3_map);
    A(k, :) = reshape(B_z, N^2, 1);
end
end
