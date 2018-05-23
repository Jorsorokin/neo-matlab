function [features,dip,dip_pval,entropy,keepCoeffs] = sparsify_wavelet_features( features,reduceDim )
    [m,nCoeffs,c] = size( features );
    dip_pval = zeros( nCoeffs,c );
    dip = dip_pval;
    entropy = dip_pval;
    for j = 1:c
        for i = 1:nCoeffs
            P = histcounts( features(:,i,j),40,'normalization','probability' );
            entropy(i,j) = -sum( P.*log( P+eps ) );
            %[dip(i,j),dip_pval(i,j)] = HartigansDipSignifTest( P,100 );
        end
    end
    
    meanEntropy = mean( entropy,2 );
    keepCoeffs = meanEntropy > quantile( meanEntropy,0.9 );
    %entropy = zscore( entropy ); % now channels are normalized
    %keepCoeffs = entropy(:) > quantile( entropy(:),0.95 );
    %keepCoeffs = dip_pval(:) < (0.05 / (nCoeffs*c));
    features = features(:,keepCoeffs,:);
    features = reshape( features,m,nnz( keepCoeffs )*c );
    %features = features(:,keepCoeffs);
    
    % perform PCA to reduce dimensionality if desired
    if reduceDim > 0 && reduceDim < size( features,2 )
        [u,s,~] = svd( features,'econ' );
        features = u(:,1:reduceDim)*s(1:reduceDim,1:reduceDim);
    end
end
            
        