function [spikeVec,trials] = spikemat2vec( spikeMat )
    % [spikeVec,trials] = spikemat2vec( spikeMat,epochLength )
    %
    % converts an N x nTrials matrix of spike times into row vector of spike times,
    % and also creates a vector "trials" consisting of numbers representing the trial (column)
    % that each spike time belongs to
    
    if size( spikeMat,1 ) == 1
        spikeVec = spikeMat;
        trials = ones( 1,numel( spikeMat ) );
        return
    end
    
    nPts = sum( ~isnan( spikeMat ) ); % # of points for each trial
    nTrials = size( spikeMat,2 );
    spikeVec = reshape( spikeMat,1,numel( spikeMat ) );
    spikeVec( isnan( spikeVec ) ) = [];
    trials = zeros( 1,numel( spikeVec ),'uint8' );
    counter = 0;
    for j = 1:nTrials
        trials(counter+1:counter+nPts(j)) = repmat( uint8(j),1,nPts(j) );
        counter = counter + nPts(j);
    end
end
    