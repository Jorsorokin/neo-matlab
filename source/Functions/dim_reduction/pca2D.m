function [V,E,proj] = pca2D( X,varargin )
% [V,E] = pca2D( X ) computes the 2-dimensional PCA of the 3D tensor X
% over the 3rd dimension of X and returns the principal components 
% (eigenvectors) and eigenvalues V and E.
%
% [V,E] = pca2D( X,nDim ) returns the top nDim principal components and
% eigenvalues. nDim must be >= 1 and less than the # of elements in the
% smallest dimension of any horizontal slice of X. 
%
% [V,E] = pca2D( X,varExp ) returns the number of principal components and
% eigenvalues that explain at least "varExp" percent of the total variance.
% "varExp" must be a value in the set (0,1). 
%
% [...,proj] = pca2D( X,... ) returns a projection of the original data in X onto the
% principal components in V
%
% X is an n x m x p tensor, with each horizontal slice X(:,:,j)
% representing the jth 2D sample (i.e. an image, multi-channel recording,
% etc.) of all p-samples. Alternatively, X may be a single sample tensor described by 
% 3 different parameters (i.e. and image over time). As noted above, the 2d
% PCA is performed over each of the jth slices of the 3rd dimension, and
% projections are computed by projecting the second dimesnion onto the
% principal components. Thus, please ensure your input tensor is of the
% correct format for your needs! 
%
% Written by Jordan Sorokin 
% jorsorokin@gmail.com
% 2/3/2018

% check inputs
[n,m,p] = size( X );

% subtract the average of X over the specified axis
Xhat = mean( X,3 );
Xbar = X - Xhat;

% create the covariance matrix of the summed slices of X
C = zeros( m,m );
for j = 1:p
    C = C + Xbar(:,:,j)'*Xbar(:,:,j); % covariance of jth sample 
end

% find the eigen vectors/values of the covariance matrix
[V,E] = eig( C );

% find # of eigen vectors to keep
switch nargin 
    case 1 
        nDim = size( V,2 );
    otherwise
        if varargin{1} >= 1
            nDim = varargin{1};
        else
            e = diag( E );
            varianceExplained = cumsum( e ) / sum( e );
            nDim = find( varianceExplained >= varargin{1},1 );
            clear e varianceExplained
        end
end
V = V(:,1:nDim);
E = E(1:nDim,1:nDim);

% compute the projection if requested
if nargout > 2 && nDim ~= m
    proj = zeros( n,nDim,p );
    for j = 1:p
        proj(:,:,j) = Xbar(:,:,j) * V;
    end
end

   
end
    