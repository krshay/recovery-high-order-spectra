% importmanopt

% Script for Figure 2(b)
rng(100);

N = 10;
Ks = [10 11 12 15 18 20 25 37 50 75 100];
NK = length(Ks);

NumIters = 120;
errs = zeros(NumIters, NK);
costs = zeros(NumIters, NK);

As = cell(NK, NumIters);
ys = cell(NK, NumIters);
xs = zeros(NumIters, N);
sig_inits = zeros(NumIters, N);

k1k2k3k4_map = calck1k2k3k4(N);

for i=1:NumIters
    display('iteration #' + string(i));
    sig = rand(N, 1) + 1j * rand(N, 1);
    x = fft(sig);
    
    xs(i, :) = x;
    T = calcTrispectrum(x, k1k2k3k4_map);
    T_flat = reshape(T, N^3, 1);
    
    for k=1:NK
        Samples = randperm(N^3, Ks(k));
        A = zeros(Ks(k) * N^3, 1);
        A(sub2ind([Ks(k), N^3], 1:Ks(k), Samples)) = 1;
        As{k, i} = reshape(A, Ks(k), N^3);
        ys{k, i} = As{k, i} * T_flat;
    end
    
    sig_init1 = rand(N, 1) + 1j * rand(N, 1);
    x_init1 = fft(sig_init1);
    
    sig_init2 = rand(N, 1) + 1j * rand(N, 1);
    x_init2 = fft(sig_init2);
    
    sig_init3 = rand(N, 1) + 1j * rand(N, 1);
    x_init3 = fft(sig_init3);
    
    for k=1:NK
        warning('off', 'manopt:getGradient:approx');
        [z1, ~, cost1] = func_optimize(x_init1, ys{k, i}, As{k, i}, k1k2k3k4_map);
        [z2, ~, cost2] = func_optimize(x_init2, ys{k, i}, As{k, i}, k1k2k3k4_map);
        [z3, ~, cost3] = func_optimize(x_init3, ys{k, i}, As{k, i}, k1k2k3k4_map);
        zs = [z1, z2, z3];
        [M, I] = min([cost1, cost2, cost3]);
        [err, shift, x_best] = calcError(x, zs(:, I));
        errs(i, k) = err * 100;
        costs(i, k) = M;
    end
end

success_rates = sum(errs < 0.005, 1) / size(errs, 1) * 100;

figure;
plot(Ks, success_rates);
xlabel('K');
ylabel('Success rate [%]');
title(['Success rate vs. K, for N = ', num2str(N), '. Estimation from samples of the trispectrum.']);
ylim([0 100]);
xlim([0 max(Ks)]);
grid on; grid minor;
