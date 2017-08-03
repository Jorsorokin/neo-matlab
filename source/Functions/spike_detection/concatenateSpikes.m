function concatSpikes = concatenateSpikes( spikes )
%
% combines the n x m x c spike matrix into an n*c x m matrix (i.e. for
% tetrode data)

[n,m,c] = size( spikes );
concatSpikes = zeros( n*c,m );

for j = 1:c
    concatSpikes([1:n] + (j-1)*n,:) = spikes(:,:,j);
end