function [S,D,neighbors] = similarityMatrix( X,K,type )
% [S,D,neighbors] = similarityMatrix( X, (K, type) )
%
% Computes the pairwise similarity matrix of the columns in X. "type" can
% be any valid distance metric as defined by "pdist.m" (refer to matlab
% documentation). "K" refers to the allowed maximum nearest neighbors to
% each point X. By default, "K" will equal the number of rows of X (all
% points)
%
% returns: 
%       S_ij = 1/K * e^( -||xi - xj||2 / 2a^2 )
%   
%   where "a" = average nearest-neighbor distance

% check inputs
[N,~] = size( X );
if nargin < 2
    K = N-1;
else
    if K >= N
        error('K must be < the number of points');
    end
end
 
if nargin < 3
    type = 'euclidean';
end

% preallocate
neighbors = zeros( N,K );
S = zeros( N,N );

% compute distance matrix S via "pdist2" function
D = pdist2( X,X,type );

% loop over each point in D, compute the K nearest neighbors of that point
for j = 1:N
    
    % get all distances to this point and sort
    [distances,idx] = sort( D(j,:),'ascend' );
    neighbors(j,:) = idx(2:K+1); 
    
    % find the average distance of nearest neighbors
    avg = mean( distances(neighbors(j,:)) ); 

    % compute S for the nearest neighbors
    S(j,neighbors(j,:)) = 1/K * exp( -bsxfun( @rdivide,distances(neighbors(j,:)),2*avg^2 ) );
end

end