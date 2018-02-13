function [labels,model] = sort_clusters( X,K,method,varargin )
% function [labels,model] = sort_clusters( X,K,method,varargin )
%
% sort the data "X" using the specified clustering method
% 
% Inputs:
%   X - n x d matrix, n = observations, d = dimensions
%
%   K - scalar (# of clusters)
%
%   method - the sorting method:
%       'EM-GMM'    <- expectation maximization guassian mixture modeling
%       'mEM-GMM'   <- masked EM-GMM
%       'EM-TMM'    <- EM t-distribution mixture modeling
%       'Km'        <- K-means
%       'VB'        <- variational bayes GMM
%       'mVB'       <- masked VB
%       'DBSCAN'    <- density based spatial clustering
%       'HDBSCAN'   <- hierarchical density based spatial clustering
%       'Spectral'  <- spectral clustering
%
%   (mask) - the c x n mask matrix, necessary for 'mEM-GMM' and 'mVB'
%
%   (sigma) - SD parameter for gaussian rbf kernel 
%             (default = 1)
%
%   (eps) - the epsilon ([0:1]) neighborhood radius (for DBSCAN or Spectral)
%           (default = 0.1)
%
%   (neighbors) - the # of neighbors for Spectral Clustering or 
%                 the k-nearest neighbor (minpts) for DBSCAN / HDBSCAN 
%                 (default = 5)
%
%   (minclustsize) - the minimum # of points in a cluster for the cluster to be kept
%                    (default = 5)
%
%   (neighborType) - type of neighborhood for spectral clustering
%                       'kNN'
%                       'mkNN'
%                       'eps'
%                   (default = 'kNN')
%
%   (kernelWeighting) - boolean flag for weighting affinity matrix by
%                       guassian RBF or leaving unweighted (euclidean)
%
%   (outlierthresh) - value between [0,1] for determining if a point is an
%                     outlier in the clustering scheme (for HDBSCAN)
%                     (default = 0.9)
%
%   * All optional arguments supplied in name-value pairs
%
% can optionally provide a mask matrix (requred for "mEM-GMM" and "mVB")

% check data size
n = size( X,1 );
if n < K
    error( 'Fewer data points than requested number of clusters' );
end

% get optional arguments
p = parse_inputs( varargin );

% cluster via different methods
fprintf( 'Clustering via %s...\n',method );
switch method
    case 'EM-GMM'
        labels = mixGaussEm( X',K );
    
    case 'mEM-GMM'   
        if any( isnan( p.mask ) )
            error( 'must supply mask matrix for maskedEM-GMM' );
        end
        labels = maskedEM_GMM( X,p.mask,K,0 ); % no search
        
    case 'EM-TMM'
        fprintf( 'EM-TMM currently under development\n. Try a different sorting routine\n' );
        labels = nan;
        
    case 'Km'
        labels = kmeans( X,K );
        
    case 'VB'
        labels = mixGaussVb( X',K );
    
    case 'mVB'
        if any( isnan( p.mask ) )
            error( 'must supply mask matrix for maskedEM-GMM' );
        end
        virtualX = compute_noise_ensemble( X,mask );
        labels = mixGaussVb( virtualX.y',K ); % uses the virtual distribution of X for clustering
    
    case 'DBSCAN'
        labels = DBSCAN( X,p.eps,p.neighbors );
        
    case 'HDBSCAN'
        hdbscan = HDBSCAN( X ); % creates an instance of the HDBSCAN cluster object
        hdbscan.run_hdbscan( p.neighbors,p.minclustsize,K,2,p.outlierThresh,false );
        labels = hdbscan.labels;
        
    case 'Spectral'
        if strcmp( p.neighborType,'eps' )
            W = adjacency_graph( X,p.neighborType,p.eps,p.kernelWeighting,p.sigma );
        else
            W = adjacency_graph( X,p.neighborType,p.neighbors,p.kernelWeighting,p.sigma );
        end
        %[labels,alpha] = wMSC( W,p.sigma,K ); 
        if p.kernelWeighting
            labels = SpectralClustering( W,K,2 );
        else
            labels = SpectralClustering( W,K,1 ); % no weighting
        end
end

% make into correct format
labels = uint8( labels );
if ~isrow( labels )
    labels = labels';
end

% remove labels (set to 0) with # pts < minclustsize
trueClusts = labels > 0;
nPts = sum( bsxfun( @minus,labels(trueClusts),repmat( unique( labels(trueClusts) )',1,nnz(trueClusts) ) )==0,2 );
smallClusters = find( nPts < p.minclustsize );
labels(ismember( labels,smallClusters )) = 0;

% get the sorting model
if nargout > 1
    switch method
        case {'EM-GMM','Km','VB'}
            model = get_sorting_model( X,labels,method );

        case 'DBSCAN'
            model = nan;
            
        case 'Spectral'
            model = nan;
%             model.W = W; % adjacency graph
%             model.D = sum( W,2 ); % degree vector
%             model.alpha = alpha; % eigen-decomp of the centered rotation: inv(D)*W
%             model.neighborType = p.neighborType; % type of neighborhood
%             model.nNeighbors = p.neighbors; % # neighbors
%             model.sigma = p.sigma;
%             model.K = K; % # of clusters
        case 'HDBSCAN'
            model = hdbscan.model;
    end
end

%% Parser
function p = parse_inputs( inputs )
    names = {'mask','sigma','eps','neighbors','neighborType',...
        'minclustsize','kernelWeighting','outlierThresh'};
    defaults = {nan,1,0.1,15,'kNN',5,0,0.9};

    p = inputParser;
    for j = 1:numel( names )
        p.addParameter( names{j},defaults{j} );
    end

    parse( p,inputs{:} );
    p = p.Results;
    
    % check for nans
    if isnan( p.eps )
        p.eps = 0.1;
    end
    
    if isnan( p.sigma )
        p.sigma = 1;
    end
    
end

end
        