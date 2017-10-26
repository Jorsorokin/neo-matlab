function [neighbors,distances,D] = compute_nearest_neighbors( X,Y,k )
% [neighbors,distances,D] = compute_nearest_neighbors( X,Y,k )
%
% comptues the nearest k neighbor distances for each point in X. 
% X = n x d matrix, with n = observations, d = variables (dimensions)
% Y = m x d matrix, with m = observations, d = variables (dimensions)
%
% neighbors, distances = n x k matrices (default to loop over X)
%
% D is a sparse n x m matrix of the pairwise distances between nearest neighbors.
% Note that it is only reasonable to output D if k << min(n,m). 

[n,~] = size( X );
[m,~] = size( Y );

% use the compliex MEX code for fast NN calculation
if any( k > [n,m] )
    error( 'requested more neighbors than number of points' );
end

% calculate the ranges for looping
neighbors = zeros( n,k );
distances = zeros( n,k );
[start,stop] = find_blockIndex_range( m,n );

for i = 1:numel( start )
    
    % euclidean distance
    pts = X(start(i):stop(i),:);
    S = compute_pairwise_dist( pts,Y );
    
    % nearest neighbors
    [S(:),idx] = sort( S,'ascend' ); % sorts smallest - largest columns for each row
    neighbors(start(i):stop(i),:) = idx(2:k+1,:)';
    distances(start(i):stop(i),:) = real( sqrt( S(2:k+1,:) ) )';
end

% compute the sparse distance matrix
if nargout > 2
    D = sparse( reshape( double( neighbors' ),n*k,1),...
                reshape( repmat( 1:m,k,1 ),m*k,1 ),...
                reshape( double( distances ),n*k,1 ),...
                n,m );
end

end