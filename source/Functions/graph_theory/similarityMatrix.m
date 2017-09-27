function S = similarityMatrix( X,K )
% S = similarityMatrix( X, (K) )
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

% preallocate
n_i = zeros( N*K,1 ); % neighbor vector for rows
n_j = n_i; % neighbor vector for columns
val = n_i; % values of neighbor similarities

% loop over each point in D, compute the K nearest neighbors of that point
for j = 1:N
    
    ind = (j-1)*K+1:j*K;
    
    % get all distances to this point and sort
    distances = pdist2( X(j,:),X );
    [~,idx] = sort( distances,'ascend' );
    neighbors = idx(2:K+1);
    n_i(ind) = j;
    n_j(ind) = neighbors;
    
    % find the average distance of nearest neighbors
    avg = mean( distances(neighbors) ); 
    
    % compute similarity
    val(ind) = 1/K * exp( -bsxfun( @rdivide,distances(neighbors),2*avg^2 ) );
end

% construct our similarity matrix
S = sparse( n_i,n_j,val,N,N );
    
end