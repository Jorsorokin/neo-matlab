function [labels,model] = sort_clusters( X,K,method,varargin )
% function [labels,model] = sort_clusters( X,K,method,(mask) )
%
% sort the projected data "X" using the specified clustering method and the
% number of clusters "K". Options for "method" include: EM-GMM, EM-TMM,
% Km, VB, DBSCAN, Spectral, and Agglomerative
%
% can optionally provide a mask matrix (requred for "mEM-GMM" and "mVB")

% check data size
if size( X,1 ) < K
    disp( 'Fewer data points than requested number of clusters' );
    labels = nan;
    return
end

% get mask if it exists
if nargin > 3 && ~isempty( varargin{1} )
    mask = varargin{1};
end

% TEMP
model = nan;
labels = nan;

if ~strcmp( method,'DBSCAN' ) && isempty( K )
    K = inputdlg( 'Number of clusters' );
    K = str2double( cell2mat( K ) );
    if isempty( K )
        return;
    end
end

% cluster via different methods
fprintf( 'Clustering via %s...\n',method );
switch method
    case 'EM-GMM'
        labels = mixGaussEm( X',K );
    
    case 'mEM-GMM'   
        labels = maskedEM_GMM( X,mask,K,0 ); % no search
        
    case 'EM-TMM'
        fprintf( 'EM-TMM currently under development\n. Try a different sorting routine\n' );
        labels = nan;
        
    case 'Km'
        % create random centers
        labels = kmeans( X,K );
        
    case 'VB'
        labels = mixGaussVb( X',K );
    
    case 'mVB'
        virtualX = compute_noise_ensemble( X,mask );
        labels = mixGaussVb( virtualX.y',K ); % uses the virtual distribution of X for clustering
    
    case 'DBSCAN'
        eps = str2double( cell2mat( inputdlg( 'Radius of neighborhood [0:1]:' ) ) );
        labels = DBSCAN( X,eps,10 );
        
    case 'Spectral'
        neighbors = inputdlg( 'Number of neighbors?' );
        neighbors = str2double( cell2mat( neighbors ) );
        if isempty( K ) || isempty( neighbors )
            labels = nan;
            return
        end
        [~,A] = adjacency_matrix( X,'kNN',neighbors );
        labels = SpectralClustering( A,K,2 );
end

% make into correct format
labels = uint8( labels );
if ~isrow( labels )
    labels = labels';
end

end
        