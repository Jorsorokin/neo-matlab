function [features,W,mapping] = compute_spike_features( X,varargin )
    % [features,W,mapping] = compute_spike_features( X,(ndim,method,mask,distances,concatenate) )
    %
    % wrapper for "compute_mapping.m" to extend the mapping method to ICA and
    % to to take masking matrices and channel maps into account
    %
    % Inputs:
    %     X: n x m x c data matrix
    %           n = # observations (spikes)
    %           m = # variables (time points) 
    %           c = # channels/electrodes
    %   
    %     (ndim): # dimensions to project onto (default = 3)
    %    
    %     (method): projection method (default = 'PCA')
    %    
    %     (mask): c x n masking matrix, with 0 <= (i,j) <= 1
    %    
    %     (distances): ch x p matrix specifying euclidean positions of each channel
    %    
    %     (concatenate): boolean indicating whether to concatenate all spikes across channels
    %                    or perform feature mapping for each channel separately
    %
    % Outputs:
    %     features: an n x ndim' matrix, where ndim' is one of
    %               (ndim, ndim*c, ndim + p, ndim*c + p), depending on the inputs
    %       
    %     W: an m x ndim (x c) projection matrix or tensor  
    %   
    %     mapping: a structure of the mapping
    %
    % by Jordan Sorokin, 8/31/17

    % check inputs
    [n,m,c] = size( X );
    mapping = struct;
    mapping.useMask = false;
    mapping.concatenate = true;
    mapping.useDistance = false;
    mapping.sizeData = [n,m,c];
    
    if nargin > 1 && ~isempty( varargin{1} )
        ndim = varargin{1};
    else
        ndim = 3;
    end
    if nargin > 2 && ~isempty( varargin{2} )
        method = varargin{2};
    else
        method = 'PCA';
    end
    if nargin > 3 && ~isempty( varargin{3} )
        mask = varargin{3};
        if all( isnan( mask(:) ) )
            clear mask
        else
            mapping.useMask = true;
        end
    end
    if nargin > 4 && ~isempty( varargin{4} ) && all( ~isnan( varargin{4}(:) ) )
        distances = varargin{4} + 1; % added 1 to avoid any distance == 0, which will affect spike-distance mapping
        mapping.useDistances = true;
    end
    if nargin > 5 && ~isempty( varargin{5} )
        concatenate = varargin{5};
        mapping.concatenate = concatenate;
    else
        concatenate = true;
    end

    % mask the channels according to the mask matrix
    if exist( 'mask','var' )
        X = maskchans( X,mask,false );
    end
    mu = mean( X,1 );
    X = X - mu;
    mapping.mu = mu;
    mapping.method = method;

    % compute regularized covariance matrix 
    if exist( 'mask','var' ) && ~concatenate && any( strcmp( method,{'PCA','pca'} ) ) 
        X2 = concatenateSpikes( X );
        mask2 = reshape( mask,numel( mask ),1 );
        reg_cov = cov( X2(mask2==0,:) );
        clear X2 mask2
    end

    % get distances to add to features
    if exist( 'distances','var' )
        if exist( 'mask','var' )
            E_spikeXY = get_spike_XY( X,distances,mask );
        else
            E_spikeXY = get_spike_XY( X,distances );
        end
    else
        E_spikeXY = [];
    end

    % PROJECTIONS
    if concatenate
        
        mapping.newDim = ndim;

        % first perform 2D PCA for dimensionality reduction, which 
        % vastly speeds up the subsequent mapping
        if size( X,3 ) > 1
            [V,E,proj] = pca2D( permute( X,[3,2,1] ),0.85 ); % 85 % variance explained
            proj = reshape( proj,c*size( proj,2 ),n )';
            mapping.pc2D.V = V;
            mapping.pc2D.E = E;
        else
            proj = X;
        end

        % compute the appropriate mapping on the reduced X
        switch method
            case {'ICA','ica'}
                [W,features,~] = fastica( proj,'numOfIC',ndim,'lasteig',ndim*2 );
                W = W';

            otherwise
                [features,map] = compute_mapping( proj,method,ndim );
                mapping.fullMap = map;
                if isfield( map,'M' )
                    W = map.M;
                else
                    W = nan;
                end
        end
    else
        
        mapping.newDim = ndim*c;
        
        W = zeros( m,ndim,c );
        features = zeros( n,ndim*c );

        % loop over channels and perform mapping for each
        for ch = 1:c
            x = X(:,:,ch);

            switch method
                case {'ICA','ica'}
                    map = fastica( x,'numOfIC',ndim,'lasteig',ndim*2 );
                    W(:,:,ch) = map';

                otherwise
                    if exist( 'mask','var' ) && any( strcmp( method,{'PCA','pca'} ) )
                        if nnz( mask(ch,:)>0 )>=2
                            covmat = cov( x ) + reg_cov;
                        else
                            covmat = reg_cov;
                        end

                        % compute PCs
                        [map,e] = eig( covmat,'vector' );
                        [~,idx] = sort( e,'descend' );
                        W(:,:,ch) = map(:,idx(1:ndim));
                        fullMap = struct;
                        fullMap.covmat = covmat;
                        fullMap.eigvec = map;
                        fullMap.eigval = e;
                        mapping.fullMap(ch) = fullMap;
                        
                    else
                        [f,mapping] = compute_mapping( x,method,ndim );
                        mapping.fullMap(ch) = mapping;
                        if isfield( mapping,'M' )
                            W(:,:,ch) = mapping.M;
                        else
                            W = nan;
                        end
                    end
            end

            % project all points from this channel onto the mapping
            thisind = ch*ndim - ndim + 1:ch*ndim;
            if exist( 'f','var' )
                features(:,thisind) = f;
            else
                features(:,thisind) = X(:,:,ch) * W(:,:,ch);
            end
        end

    end

    % add distances
    if isfield( mapping,'conn_comp' ) && ~isempty( E_spikeXY )
        E_spikeXY = E_spikeXY(mapping.conn_comp,:);
    end
    features = [features,E_spikeXY];
end
