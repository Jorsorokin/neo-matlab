function [Cinv,L] = compute_noiseCov_inv( C,alpha,nPts )
    % [Cinv,L] = compute_noiseCov_inv( C,alpha,nPts )
    %
    % Compute the inverse of the noise covariance matrix C efficiently using 
    % cholesky decomposition. "alpha" controls the relative contribution of
    % cross-electrode covariance via the equation:
    %
    %       alpha*C + (1-alpha)*diag(C)
    %
    
    % re-weight the noise covariance
    n = size( C,1 );
    nChan = n / nPts;
    diagC = zeros( n,n );
    startPt = 1;
    for j = 1:nChan
        inds = startPt:startPt + nPts-1;
        diagC(inds,inds) = C(inds,inds);
        startPt = startPt + nPts;
    end

    C = alpha*C + (1-alpha)*diagC;
    
    % find inverse via cholesky decomp
    [L,p] = chol( C,'lower' );
    if p ~= 0
        fprintf( 'C is poorly conditioned!\n' );
        Cinv = nan;
        return
    end
    
    u = L \ eye( size( L,1 ) );
    Cinv = L' \ u;
end
    