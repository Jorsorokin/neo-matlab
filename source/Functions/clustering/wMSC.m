function [label,alpha] = wMSC( W,sigma,k,varargin )
% label = wMCS( W,D,sigma,k,(alpha) )
%
% Performs weighted kernel-PCA multisignal spectral clustering
%
% Inputs:
%	W - N x N undirected (symmetric) graph matrix
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
	error( 'Kernel matrix must be symmetric' );
end

% create indicator and inv( D ) variables
indicator = ones( N,1 );
d = sum( W,2 );
d(d==0) = eps;
invD = spdiags( 1./(d.^0.5),0,N,N );

% compute alpha if needed
if ~exist( 'alpha','var' )
	M = eye( N ) - (1 / (indicator'*invD*indicator)) * (indicator*indicator'*invD); % centering matrix
	if N <= 2000
		[alpha,lambda] = eig( invD*M*W,'vector' );
	else
		[alpha,lambda] = eigs( invD*M*W,k );
        lambda = diag( lambda );
	end

	% sort alpha
	[~,ind] = sort( lambda,'descend' );
	alpha = alpha(:,ind(1:k-1)); % keep only k-1 top eigenvectors
end

% compute the bias term
b = -(1 / (indicator'*invD*indicator))*(indicator'*invD*W*alpha);

% compute the scores (projections)
z = bsxfun( @plus,W*alpha,b );

% binarize the scores and eigenvectors
z_binary = sign( z ) > 0;
alpha_binary = sign( alpha ) > 0;

% create the codebook based on k-most occurring patterns of alpha_binary
[codes,codeIDX] = unique( alpha_binary,'rows' );
codeCount = get_hamming( alpha_binary,codes );
codeCount = sum( codeCount == 0 );
[~,sortedCodeIDX] = sort( codeCount,'descend' );
codebook = alpha_binary(codeIDX(sortedCodeIDX(1:k)),:);

% compute the hamming distance of each codebook vector with the binary scores
distances = get_hamming( z_binary,codebook );

% assign cluster labels to smallest distances
[~,label] = min( distances,[],2 ); 

end

%%
function dist = get_hamming( A,B )
    [n,~] = size( A );
    p = size( B,1 );
    dist = zeros( n,p );

    for j = 1:p
        dist(:,j) = sum( A ~= B(j,:),2 );
    end
end
