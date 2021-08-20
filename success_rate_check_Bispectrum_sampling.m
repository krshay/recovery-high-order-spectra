% Script for Figure 2(a)
rng(1);

N = 30;
Ks = [30 31 32 35 40 50 62 75 100 125 150];
NK = length(Ks);

NumIters = 100;
errs = zeros(NumIters, NK);
costs = zeros(NumIters, NK);

As = cell(NK, NumIters);
ys = cell(NK, NumIters);
xs = zeros(NumIters, N);

k1k2k3_map = calck1k2k3(N);

parfor i=1:NumIters
    display('iteration #' + string(i));
    sig = randn(N, 1);
    x = fft(sig);
    
    xs(i, :) = x;
    B = calcBispectrum(x, k1k2k3_map);
    B_flat = reshape(B, N^2, 1);
    
    for k=1:NK
        Samples = randperm(N^2, Ks(k));
        A = zeros(Ks(k) * N^2, 1);
        A(sub2ind([Ks(k), N^2], 1:Ks(k), Samples)) = 1;
        As{k, i} = reshape(A, Ks(k), N^2);
        ys{k, i} = As{k, i} * B_flat;
    end
    
    sig_init1 = randn(N, 1);
    x_init1 = fft(sig_init1);
    
    sig_init2 = randn(N, 1);
    x_init2 = fft(sig_init2);
    
    sig_init3 = randn(N, 1);
    x_init3 = fft(sig_init3);
    
    for k=1:NK
        [z1, cost1] = optimize([real(x_init1); imag(x_init1)], ys{k, i}, As{k, i}, k1k2k3_map);
        [z2, cost2] = optimize([real(x_init2); imag(x_init2)], ys{k, i}, As{k, i}, k1k2k3_map);
        [z3, cost3] = optimize([real(x_init3); imag(x_init3)], ys{k, i}, As{k, i}, k1k2k3_map);
        zs = [z1, z2, z3];
        [M, I] = min([cost1, cost2, cost3]);
        Z = zs(:, I);
        z = Z(1:N) + 1j * Z(N+1:end);
        err = calcError(x, z);
        errs(i, k) = err * 100;
        costs(i, k) = M;
    end
end

success_rates = sum(errs < 0.005, 1) / size(errs, 1) * 100;

figure;
plot(Ks, success_rates);
xlabel('K');
ylabel('Success rate [%]');
ylim([0 100]);
xlim([0 max(Ks)]);
grid on; grid minor;

save data_2a.mat
