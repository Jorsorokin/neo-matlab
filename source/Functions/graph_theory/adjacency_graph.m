function [D,A,L] = adjacency_graph( X,type,val,varargin )
% function [D,A,L] = adjacency_graph( X,type,val,(weight,sigma) );
% 
% computes the similarity (distance), affinity, and normalized laplacian graphs
% from the data in X. Implements batch looping and sparse matrices where
% applicable for fast and memory efficient computations. Actual similarity 
% matrices are computed using compiled MEX files for further speedups.
%
% Inputs:
%   X           - an n x m matrix, with n = observations, m = variables (dimensions)
%
%   type        - a string specifying the type of graph to compute. Valid options are:
%                   'kNN'   -  k-nearest neighbors
%                   'mkNN'  -  mutual k-nearest neighbors
%                   'eps'   -  epsilon-radius neighborhood
%
%   val         - an integer or decimal specifying:
%                   1) # of neighbors (for 'kNN' and 'mkNN')
%                   2) radius of neighborhood ball (for 'eps')
%
%   (weight)    - a boolean (true/false) indicating whether to weight the graph or not 
%                    default = false
%
%   (sigma)     - number specifying gaussian variance for weighting 
%                    default = 1
%
% Outputs:
%   D   - the similarity (distance) graph:
%
%                            | euclidean( x_i,x_j ): x_i is neighbor of x_j
%                   D(i,j) = |
%                            | 0: otherwise
%
%   A   - the affinity graph:
%
%                            | K( D(i,j) ): x_i is neighbor of x_j, and K is an RBF kernel (if "weight" is true)
%                   A(i,j) = |
%                            | 0: otherwise
%   
%   
%   L   - the normalized 'random walk' laplacian:
%
%                   L = inv(W) * (W - A), where w_i = diag{ SUM_j{ A_i } }, and W := diag( w )       
%
% Written by Jordan Sorokin, 9/16/2017

% check inputs
if nargin > 3 && ~isempty( varargin{1} )
    weight = varargin{1};
else
    weight = false;
end
if nargin > 4 && ~isempty( varargin{2} )
    sigma = varargin{2};
else
    sigma = 1;
end
[n,~] = size( X );

% compute the (sparse) graph matrix D:
switch type
    case {'kNN','mkNN'}
        if ~isnumeric( val )
            error( 'must input an integer value for k-nearest neighbor graphs' );
        end

        [~,~,D] = compute_nearest_neighbors( X,X,val );

        if strcmp( type,'kNN' )
            D = max( D,D' );
        else
            D = min( D,D' );
        end

    case 'eps'
        if (val > 1) || (val < 0)
            error( 'epsilon radius must be between [0,1]' );
        end

        D = compute_epsilon_neighborhood( X,X,val );
        D = max( D,D' ); % makes it symmetric

    otherwise
        error( 'unrecognized graph type. Please see documentation' );
end

% get the affinity matrix
if nargout > 1
    if (weight==true) && (sigma > 0)
        gaussFunc = @(X,sigma)(exp( -X.^2 ./ (2*sigma^2) )); % RBF kernel
        A = spfun( @(D)gaussFunc( D,sigma ),D );
    else
        A = spfun( @(D)gt( D,0 ),D );
    end
end

% compute the nromalized laplacian
if nargout > 2
    w = sum( A,2 );
    W = sparse( 1:n,1:n,w );
    L = W - A; % un-normalized laplacian
    w(w==0) = eps;
    W = spdiags( 1./w,0,n,n ); % inv( W )
    L = W * L; % random-walk normalization
end

end
