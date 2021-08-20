function [f] = cost_Tri(z, y, A, k1k2k3k4_map)
%cost_Tri Calculates the cost and grad for the optimization
%using the Trispectrum

%   z is the currentf guess for the fft of the signal
%   y is the measurement
%   A is the sensing matrix
%   k1k2k3k4_map is a N^3*4 matrix containing the quadrupels of frequencies

N = size(z, 1) / 2;

T = calcTrispectrum(z(1:N) + 1j * z(N+1:end), k1k2k3k4_map);

f = (1 / 2) * sum(abs(y - A * T(:)).^2);
end