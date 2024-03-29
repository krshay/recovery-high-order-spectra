function A = calc_mat_trispectrum(K, N, k1k2k3k4_map)
%calc_mat_trispectrum Calculates the sensing matrix that consists
%of trispectra, for experiment 3b

%   K is the number of samples
%   N is the length of the signal
%   k1k2k3_map is a N^3*4 matrix containing the quadruples of frequencies
A = zeros(K, N^3);
for k=1:K
    sig = 1 * (randn(N, 1) + 1j * randn(N, 1));
    z = fft(sig);
    T_z = calcTrispectrum(z, k1k2k3k4_map);
    A(k, :) = reshape(T_z, N^3, 1);
end
end
