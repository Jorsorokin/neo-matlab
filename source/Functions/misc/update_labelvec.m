function newLabels = update_labelvec( newLabels,oldLabels,maxID )
    % newLabels = update_IDvec( newLabels,oldLabels,maxID )
    %
    % given a particular label vector for data points belonging to a
    % previous cluster, updates the label vector to match the current max #
    % ID in the dataset and the previous ID of cluster that has been split
    zeroLabels = newLabels == 0;
    sameLabels = ismember( newLabels,unique( oldLabels ) );
    changeLabels = ~sameLabels & ~zeroLabels;
    newLabels(changeLabels) = newLabels(changeLabels) - min( newLabels(changeLabels) ) + maxID + 1;
end