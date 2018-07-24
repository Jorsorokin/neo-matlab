function [xcg,lags] = correlogram( train1,train2,binwidth,maxlag )
    % [xcg,lags] = correlogram( train1,train2,binwidth,maxlag )
    %
    % compute the cross-correlogram between two spike trains. If train1 and train2 are matricies,
    % the must be column-major (trials = columns, time points = rows), and have the same number of
    % trails (columns).  
    
    [nPts,nTrials] = size( train1 );
    lags = -maxlag:binwidth:maxlag;
    
    if nTrials > 1 && nPts > 1
        nBins = numel( lags ) - 1;
        xcg = zeros( 1,nBins );
        for j = 1:nTrials
            dT = train1(:,j) - train2(:,j)';
            dT(dT==0) = nan;
            xcg = xcg + histcounts( dT,lags,'Normalization','pdf' );
        end
    else
        dT = train1 - train2';
        dT(dT==0) = nan;
        xcg = histcounts( dT,lags,'Normalization','pdf' );
    end
    
    xcg = xcg * (1/nTrials);
    lags = 0.5*(lags(1:end-1)+lags(2:end));
end
    