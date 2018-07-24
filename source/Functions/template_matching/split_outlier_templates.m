function newLabels = split_outlier_templates( templates,labels,spikes,mask,nDIM,maxK )
    % newLabels = split_outlier_templates( templates,labels,spikes,mask,nDIM,maxK,minClustSize )
    % 
    % checks for outlier templates (see check_template_outliers) then
    % splits by finding the optimum # of clusters via iterative em-gmm (see
    % fit_optimum_gmm).
    %
    % newLabels will contain the updated labels from the splits 
    
    outliers = check_template_outliers( templates );
    IDX = ismember( labels,outliers );
    newLabels = resort_spikes( spikes(:,IDX,:),mask(:,IDX),labels(IDX) );
    maxID = max( labels(:,1) );
    
    for k = outliers
        idx = labels==k;
        if nnz( idx ) < 10
            continue
        end

        proj = compute_spike_features( permute( spikes(:,idx,:),[2,1,3] ),nDIM,...
                                       'PCA',mask(:,idx),[],true );
                                   
        [~,~,newLabels] = fit_optimum_gmm( proj,maxK,...
                                           'replicates',5,...
                                           'regularization',0.01,...
                                           'outlierpdf',1e-6 );                           
        labels(idx) = update_labelvec( single( newLabels ),k,maxID );
        maxID = max( labels );
    end
end