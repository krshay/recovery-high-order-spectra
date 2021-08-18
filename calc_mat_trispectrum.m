function A = calc_mat_trispectrum(K, N, k1k2k3k4_map)

A = zeros(K, N^3);
for k=1:K
    sig = 1 * (randn(N, 1) + 1j * randn(N, 1));
    z = fft(sig);
    T_z = calcTrispectrum(z, k1k2k3k4_map);
    A(k, :) = reshape(T_z, N^3, 1);
end
end
