function C = compute_noise_covariance( data,spikeTimes,totalPts,nIter )
    % C = compute_noise_covariance( data,spikeTimes,totalPts,nIter )
    %
    % Computes the block-covariance of noise between all pairs of
    % electrodes in the N x M matrix "data"
    
    [n,m] = size( data );
    msquared = m^2; 
    bigSamples = totalPts*m;
    C = zeros( bigSamples,bigSamples );
    
    for j = 1:nIter
        startIDX = randperm( n,1 );
        stopIDX = startIDX + totalPts - 1;
        while any( ismember( spikeTimes,startIDX:stopIDX ) )
            startIDX = randperm( n,1 );
            stopIDX = startIDX + totalPts - 1;
        end  
        noiseSeg = data(startIDX:stopIDX,:);
        temp = xcorr( noiseSeg ); % totalSamples*2 - 1 x msquared cross-correlation matrix
        rowCounter = 1;
        colCounter = 1;
        for chan = 1:msquared
            T = toeplitz( temp(totalPts:end,chan) );
            rowInds = rowCounter:rowCounter+totalPts-1;
            colInds = colCounter:colCounter+totalPts-1;
            C(rowInds,colInds) = C(rowInds,colInds) + T;
            bottomRow = (rowCounter+totalPts < bigSamples);
            rowCounter = rowCounter*bottomRow + max( totalPts*bottomRow,1 );
            colCounter = colCounter + (totalPts*~bottomRow);
        end
    end
    
    C = C * (1/nIter);
end