function [err, shift, x_best] = calcError(x, x_est)
%calcError Calculates the relative error using eq. (4.3)

%   x is a complex N*1 vector
%   x_est is a complex N*1 vector; the estimated x from the optimization
%   procedure

N = size(x, 1);

errs = zeros(N, 1);
x_shifted = zeros(N, N);

% shifting to account for the circular shift symmetry
for shift=0:N-1
    x_shifted(shift+1, :) = x_est.' .* exp(1j * (0:N-1) * (2 * pi() / N) * shift);
    errs(shift+1) = norm(x - x_shifted(shift+1, :).', 2) / norm(x, 2);
end

[err, idx] = min(errs);
shift = idx - 1;
x_best = x_shifted(idx, :).';
end