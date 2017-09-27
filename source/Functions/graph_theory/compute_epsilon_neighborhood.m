function D = compute_epsilon_neighborhood( X,Y,eps )
% D = compute_epsilon_neighborhood( X,Y,eps )
%
% comptues the neighbors for each point in X based on the radius of 'eps'
% X = n x d matrix, with n = observations, d = variables (dimensions)
% Y = m x d matrix, with m = observations, d = variables (dimensions)
%
% D is an n x m sparse matrix of the distances between points within 
% an 'eps' radius of one another.

[n,~] = size( X );
[m,~] = size( Y );

% calculate the ranges for looping
[start,stop] = find_blockIndex_range( m,n );
D = zeros( n,m );

for i = 1:numel( start )
    
    % euclidean distance
    pts = X(start(i):stop(i),:);
    S = sqrt( abs( bsxfun( @plus,dot( pts,Y,2 ),...
                            bsxfun( @minus,dot( pts,Y,2 ),2*(pts*Y') ) ) ) ); 

    % add to the "D" matrix
    D(start(i):stop(i),:) = S .* (S<=eps);
end

D = sparse( D );

end