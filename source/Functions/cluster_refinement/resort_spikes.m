function labels = resort_spikes( spikes,mask,labels,maxK,nDim,outlierPDF )
    % labels = resort_spike( spikes,mask,labels,maxK,nDim,outlierPDF)
    %
    % resorts the clusters of spikes (defined by the 'labels' vector) using
    % iterative EM-GMM to find the optimal number of clusters. 
    %
    % WARNING: Each cluster in the labels vector will be resorted!
    
    uID = unique( labels );
    if ~isrow( uID )
        uID = uID';
    end
    
    maxID = max( labels );
    for k = uID
        idx = labels==k;
        if nnz( idx ) < nDim
            continue
        end

        proj = compute_spike_features( permute( spikes(:,idx,:),[2,1,3] ),nDim,...
                                       'PCA',mask(:,idx),[],true );
                                   
        [~,~,newLabels] = fit_optimum_gmm( proj,maxK,...
                                           'replicates',5,...
                                           'regularization',0.01,...
                                           'outlierpdf',outlierPDF );  
                                                   
        labels(idx) = update_labelvec( single( newLabels ),labels(idx),maxID );
        maxID = max( labels );
    end
end
    