function [mu,sigma] = get_cluster_description( X,labels )
    % function [mu,sigma] = get_cluster_description( X,labels )
    %
    % returns the mean (mu) and covariances (sigma) of data in X for each
    % cluster defined by "labels"
    uID = unique( labels );
    if ~isrow( uID )
        uID = uID';
    end
    nID = numel( uID );
    nDim = size( X,2 );
    sigma = cell(1,nID);
    mu = zeros( nID,nDim );
    for id = 1:numel( uID )
        thisID = labels==uID(id);
        mu(id,:) = mean( X(thisID,:) );
        sigma{id} = cov( X(thisID,:) );
    end
end