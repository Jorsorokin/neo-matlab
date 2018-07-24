function labels = remove_small_clusters( labels,minSize )
    % labels = remove_small_clusters( labels,minSize )
    %
    % sets labels = 0 if the number of points for that cluster are less
    % than minSize
    
    uID = unique( labels(labels>0)' );
    for k = uID
        idx = (labels(:,1) == k);
        if nnz( idx ) < minSize
            labels(idx) = 0;
        end
    end
end