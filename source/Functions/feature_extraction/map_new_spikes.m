function features = map_new_spikes( spikes,W,mapping,varargin )
    % features = map_new_spikes( spikes,W,mapping,(mask,chanDistances) ) transforms new spikes
    % via the linear transformation matrix W
    %
    % Inputs:
    %   spikes - an n x m x c matrix, with:
    %       n = # spikes
    %       m = # points / spike
    %       c = # of channels
    % 
    %   W - the transformation matrix
    %
    %   mapping - structure computed from "compute_spike_features.m"
    %
    %   (mask) - a c x n masking matrix, to mask channels if
    %            "mapping.useMask" is true
    %
    %   (chanDistances) - a c x p matrix of channel locations to add to the
    %                     feature matrix if "mapping.useDistance" is true
    %
    % Outputs:
    %   features - an n x k matrix, with k = "mapping.newDim"
    %
    % Written by Jordan Sorokin, 5/16/18
    
    % check inputs
    [n,m,c] = size( spikes );
    assert( m == mapping.sizeData(2),'original dimensionality does not match mapping' );
    assert( c == mapping.sizeData(3),'original # of channels does not match mapping' );
    if nargin > 3 && ~isempty( varargin{1} )
        mask = varargin{1};
        assert( all( size( mask ) == [c,n] ),'mask must be of size c x n, and equal data size' );
    elseif mapping.useMask
            error( 'original mapping used a masking matrix, but none provided here!' );
    end
    if nargin > 4 && ~isempty( varargin{1} )
        distances = varargin{1};
        assert( size( distances,1 ) == c,'# of rows of the distance matrix must = c' );
    elseif mapping.useDistance
        warning( 'original feature matrix contained spike XY location estimation, but no distance matrix provided here!' );
    end
    
    % compute distances if distance provided
    if exist( 'distances','var' ) && mapping.useDistance
        if mapping.useMask
            E_spikexy = get_spike_xy( spikes,distances,mask );
        else
            E_spikexy = get_spike_xy( spikes,distances );
        end
    else
        E_spikexy = [];
    end
    
    % mask the data if mask provided
    if mapping.useMask
        spikes = maskchans( spikes,mask );
    end
    
    % transform the data 
    spikes = spikes - mapping.mu; % remove mean of original dataset
    if mapping.concatenate
        [~,~,proj] = pca2D( permute( spikes,[2,3,1] ),mapping.pc2D.V );
        features = proj' * W;
    else
        if ~isempty( which( 'mmx' ) )
            features = reshape( mmx( 'mult',double( spikes ),double( W ) ),n,mapping.newDim );
        else
            features = zeros( n,mapping.newDim / c,c );
            for ch = 1:c
                features(:,:,ch) = spikes(:,:,ch) * W(:,:,ch);
            end
            features = reshape( features,n,mapping.newDim );
        end
    end
    
    % add distances
    features = [features,E_spikexy];
end

    