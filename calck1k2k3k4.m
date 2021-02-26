function k1k2k3k4_map = calck1k2k3k4(N)
%calck1k2k3k4 Calculates the quadruple of frequencies, utilized in
%calcTrispectrum

%   N is the length of the signal
kmap = 0:N-1;
[X, Y, Z] = meshgrid(kmap, kmap, kmap);
k3map = [X(:) Y(:) Z(:)];
k1k2k3k4_map = zeros(size(k3map, 1), 4);
for i=1:size(k3map, 1)
    k1 = k3map(i, 1);
    k2 = k3map(i, 2);
    k3 = k3map(i, 3);
    k4 = mod(-k1 - k2 - k3, N);
    k1k2k3k4_map(i, :) = [k1 k2 k3 k4];
end
end

