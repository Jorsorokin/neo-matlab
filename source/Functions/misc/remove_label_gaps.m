function labels = remove_label_gaps( labels )
    % labels = remove_label_gaps( labels )
    %
    % removes jumps in labels (i.e. [1 1 5 5 5] becomes [1 1 2 2 2])
    
    uID = unique( labels(labels>0) );
    fixLabels = true( 1,numel( labels ) );
    if ~isrow( uID )
        uID = uID';
        fixLabels = fixLabels';
    end
    
    nID = numel( uID );
    for i = 1:nID
        idx = (labels == uID(i) & fixLabels);
        labels(idx) = i;
        fixLabels(idx) = false;
    end
end
        
        
        