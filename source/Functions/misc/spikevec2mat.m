function spikeMat = spikevec2mat( spikeVec,trials )
    % spikeMat = spikevec2mat( spikeVec,trials )
    %
    % converts a spike-time vector with spiketimes from multiple trials into an N x nTrials matrix, 
    % where N represents the maximum # of spike times in any trial. Trials with fewer spike times
    % than N will be padded with NaN. spikeVec and trials should both be either row or column vectors
    
    n = numel( trials );
    nTrials = max( trials );
    if isrow( trials )
        trials = trials';
        spikeVec = spikeVec';
    end

    trialIDX = (trials == repmat( 1:nTrials,n,1 ));
    nPts = sum( trialIDX );
    spikeMat = nan( max( nPts ),nTrials );
    
    for trial = 1:nTrials
        spikeMat(1:nPts(trial),trial) = spikeVec(trialIDX(:,trial));
    end
end