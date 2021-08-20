function A = calc_mat_bispectrum(K, N, k1k2k3_map)

A = zeros(K, N^2);
for k=1:K
    sig = randn(N, 1) + 1j * randn(N, 1);
    z = fft(sig);
    B_z = calcBispectrum(z, k1k2k3_map);
    A(k, :) = reshape(B_z, N^2, 1);
end
end
