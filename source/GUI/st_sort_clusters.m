function [labels,model,probabilities] = st_sort_clusters( X,K,method )
% function [labels,model,probabilities] = st_sort_clusters( X,K,method )
%
% sort the projected data "X" using the specified clustering method and the
% number of clusters "K". Options for "method" include: EM-GMM, EM-TMM,
% Km, VB

% check data size
if size( X,1 ) < K
    disp( 'Fewer data points than requested number of clusters' );
    labels = nan;
    model = nan;
    return;
end

% cluster via different methods
fprintf( 'Clustering via %s...\n',method );
model = struct;
probabilities = nan;
switch method
    case 'EM-GMM'
        [labels,model,~,probabilities] = mixGaussEm( X',K );
    case 'EM-TMM'
        % do things 
    case 'Km'
        [labels,mu] = kmeans( X,K );
        model.mu = mu;
        [~,model.sigma] = get_cluster_description( X,labels );
    case 'VB'
        % do things
end

% make into correct format
labels = uint8( labels );
if ~isrow( labels )
    labels = labels';
end

end
        