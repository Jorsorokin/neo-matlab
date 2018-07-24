function spikes = expandSpikes( spikes,nChan )
    % spikes = expandSpikes( spikes,nChan )
    %
    % re-expands the n*c x m matrix of spike waveforms into an n x m x c
    % tensor (undoing the results of "concatenateSpikes" )
    
    [N,m] = size( spikes );
    n = N / nChan;
    spikes = permute( reshape( spikes,n,nChan,m ),[1,3,2] );
end