function label = wMSC( W,D,sigma,k,varargin )
% label = wMCS( W,D,sigma,k,(alpha) )
%
% Performs weighted kernel-PCA multisignal spectral clustering
%
% Inputs:
%	W - N x N undirected (symmetric) graph matrix
%	
%	D - N x N diagonal degree matrix 
%	
%	sigma - the SD of the gaussian rbf kernel
% 	
%	k - the # of clusters
%
% 	(alpha) - the weighting (projection) vector. Providing 'alpha'
%			  allows for clustering of new points not originally used 
%			  in the eigen-decomposition of the kernel matrix
%
% Outputs:
%	label - 1 x N vector of cluster labels for each point 1:N
%
% By Jordan Sorokin, 9/26/2017

if nargin > 4 && ~isempty( varargin{1} )
	alpha = varargin{1};
end

% check size of the kernel matrix and degree matrix
[N,M] = size( W );
if N ~= M
	error( 'Kernel matrix W must be symmetric' );
end
if N ~= size( D,1 )
	error( 'Kernel and degree matrices must have equal sizes' );
end

% compute alpha if needed
if ~exist( 'alpha','var' )
	M = eye( N ) - (1 / (indicator'*invD*indicator))*(indicator*indicator'*invD); % centering matrix
	if N <= 2000
		[alpha,lambda] = eig( invD*M*W );
	else
		[alpha,lambda] = eigs( invD*M*W );
	end

	% sort alpha
	[~,ind] = sort( diag( lambda ),'descend' );
	alpha = alpha(:,ind(1:k-1)); % keep only k-1 top eigenvectors
end

% compute the bias term
b = -(1 / (indicator'*invD*indicator))*(indicator'*invD*W*alpha);

% compute the scores (projections)
z = bsxfun( @plus,W*alpha,b );

% binarize the scores and eigenvectors
z_binary = sign( z );
alpha_binary = sign( alpha );

% create the codebook based on k-most occurring patterns of alpha_binary
[~,codeIDX] = unique( alpha_binary','rows' );
codes = alpha_binary(codeIDX,:);
nCodes = numel( codeIDX );
codeCount = zeros( 1,nCodes );
codeCount = pdist2( alpha_binary,codes,'hamming' );
codeCount = sum( codeCount==0 ); 
% for j = 1:nCodes
% 	codeCount(j) = nnz( ismember( alpha_binary,alpha_binary(codeIDX(j),:) ) );
% end
[~,sortedCodeIDX] = sort( codeCount,'descend' );
codebook = alpha_binary(codeIDX(sortedCodeIDX(1:k)),:);

% compute the hamming distance of each codebook vector with the binary scores
distances = pdist2( z_binary,codebook,'hamming' );

% assign cluster labels to smallest distances
[~,labels] = min( distances,[],2 ); 

end
