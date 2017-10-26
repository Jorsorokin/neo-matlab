function [neighbors,distances] = nn_batch( X,Y,k )

[n,~] = size( X );
[m,~] = size( Y );
if any( k > [n,m] )
    error( 'requested more neighbors than number of points' );
end

% calculate the ranges for looping
neighbors = zeros( n,k );
distances = zeros(n,k );
[start,stop] = find_blockIndex_range( m,n );

for i = 1:numel( start )
    
    % euclidean distance
    pts = X(start(i):stop(i),:);
    S = compute_pairwise_dist( pts,Y );
    
    % nearest neighbors
    [S(:),idx] = sort( S,'ascend' ); % sorts smallest - largest columns for each row
    neighbors(start(i):stop(i),:) = idx(2:k+1,:)';
    distances(start(i):stop(i),:) = S(2:k+1,:)';
end
end