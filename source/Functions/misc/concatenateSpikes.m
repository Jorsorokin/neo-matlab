function concatSpikes = concatenateSpikes( spikes )
%
% combines the n x m x c spike matrix into an n*c x m matrix (i.e. for
% tetrode data). n = points, m = observations, c = channels

[n,m,c] = size( spikes );
concatSpikes = reshape( permute( spikes,[1,3,2] ),n*c,m );

end