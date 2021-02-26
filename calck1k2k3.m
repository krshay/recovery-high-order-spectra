function k1k2k3_map = calck1k2k3(N)
%calck1k2k3 Calculates the triplets of frequencies, utilized in
%calcBispectrum

%   N is the length of the signal
kmap = 0:N-1;
[X, Y] = meshgrid(kmap, kmap);
k2map = [X(:) Y(:)];
k1k2k3_map = zeros(size(k2map, 1), 3);
for i=1:size(k2map, 1)
    k1 = k2map(i, 1);
    k2 = k2map(i, 2);
    k3 = mod(-k1 - k2, N);
    k1k2k3_map(i, :) = [k1 k2 k3];
end
end

