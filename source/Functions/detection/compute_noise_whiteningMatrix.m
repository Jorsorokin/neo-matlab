function [C,W] = compute_noise_whiteningMatrix( spikes,mask,noiseThresh )
    % [C,W] = compute_noise_whiteningMatrix( spikes,mask,noiseThresh )
    %
    % Compute the noise whitening matrix across all i,j electrode pairs
    % using the eigen decomposition of the i,j noise covariance matrix.
    %
    % THe noise covariance is computed from segments of data in which 
    % a given spike k is not present on either electrode i or j, and in
    % addition there is not a large deflection on either electrode for
    % that particular snippet of data, which is controlled via the
    % "noiseThresh" input
    %
    % Inputs:
    %   spikes - an n x m x c tensor, with n = pts, m = spikes, c = chans
    %
    %   mask - an m x c masking matrix
    %
    %   noiseThresh - a positive scalar indicating the threshold to
    %                 consider a particular snippet of data contaminated
    %                 with a signal other than noise
    %
    % Outputs:
    %   C - an n*c x n*c covariance matrix of inter-electrode covariances
    %
    %   W - an n*c x n*c whitening matrix computed from C
    %
    % By Jordan Sorokin, 5/30/18
    
    [n,~,c] = size( spikes );
    totalSamples = n*c;
    W = zeros( totalSamples,totalSamples );
    C = zeros( totalSamples,totalSamples );
    rowCounter = 1;
    colCounter = 1;

    for i = 1:c
        idx = (mask(i,:) == 0) & (~any( abs( spikes(:,:,i) ) > noiseThresh ));
        colInds = colCounter:colCounter + n - 1;

        for j = 1:c
            rowInds = rowCounter:rowCounter + n - 1;
            
            idx2 = idx & (mask(j,:) == 0) & (~any( abs( spikes(:,:,j) ) > noiseThresh )); 
            spikes_i = spikes(:,idx2,i);
            spikes_j = spikes(:,idx2,j);

            C_ij = cov( [spikes_i,spikes_j]' ) + eye( n )*eps;
            [u,s] = eig( C_ij );
            W_ij = u*diag(1./(diag(s) + eps))*u';

            W(rowInds,colInds) = W_ij;
            C(rowInds,colInds) = C_ij;
            
            bottomRow = (rowCounter + n < totalSamples);
            rowCounter = rowCounter*bottomRow + max( n * bottomRow,1 );
        end
        
        colCounter = colCounter + n;
    end
end
