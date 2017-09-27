function [ID,Params] = spikesort( spikes,Fs,varargin )
% ------------ [ID,Params] = spikesort( spikes,Fs,varargin ) ---------
%
%   Sort columns of extracellularly recorded spikes by first projecting
%   spikes onto a smaller subspace for extracting informative features of
%   the spike waveform shapes (see "method" below), then using a specified 
%   clustering algorithm (see "sortMethod" below) to cluster into groups
%
%   This function can also work in a parallel computing environment
%
%                   >>> INPUTS >>>
%
%   spikes: m x n x c matrix of spike snips with m = length of each spike waveform,
%               n = # of spikes, c = # of channels
%
%   Fs: the sampling rate (used for upsampling and width calculation)
%           * All optional inputs should be given as name-value pairs *
%
%   (method): string specifying the feature-selection method to use for
%           sorting. Valid arguments are "raw", "ICA", and any valid decomposition
%           technique outlined in "compute_mapping.m"
%
%   (decomp): true/false...if true, will perform dimensionality reduction
%           via kPCA on the extract waveform features (amplitude, ahp, etc). 
%           Note, this only applies if "raw" is chosen as "method"  
%   
%   (init): the number of clusters for EM algorithm or a model structure from
%         a previous run of an EM algorithm (containing the mu, sigma, and 
%         weights for each cluster). Default is to set to 2 clusters and
%         not use a previous model.
%   
%   (mapping): an n x m projection matrix to be mapped onto.
%
%   (features): an n x m matrix of features computed from a previous spike
%             sorting run. If provided, no PCs will be computed
%
%   (level): the total number of featuresto keep for clustering. 
%          If "level" is specified, the algorithm will keep the top "level" number of features.
%          (default = 5)
%
%   (reject): the min probability of point X belonging to any cluster.
%           probabilities lower than this will result in rejection of point
%           X from the cluster analysis (i.e. ID = 0). Default value is 0.6
%           (can be any value in the set [0,1.0]) 
%   
%   (plotting): true/false...plots the sorted spike snips and scatter plot of
%             first three features or PC space (default = 0).
%             Additionally plots the within/between distance and R value if
%             "search" is set to 1. 
%   
%   (par): true/false...if 1, sets the "UseParallel" option to true for Kmeans
%          and invokes parallel computing environment
%   
%   (search): true/false...if 1, searches for an optimal cluster number using a
%           GMM distribution
%
%   (mask): a sparse or full c x n masking matrix (see "double_flood_fill.m")
%
%   (sortMethod): a string specifying the sorting method to use. Valid options are
%                   "EM-GMM" (expectation maximization for gaussian mixture models) <-- default
%                   "mEM-GMM" (masked EM-GMM)
%                   "VB-GMM" (variational bayes for gaussian mixture model)
%   
%   (distances): a vector or matrix of channel distances for multi-channel spike matrices
%         
%                   <<< OUTPUTS <<<
%
% ID: n x 1 vector of cluster ID's 
%
% Params: strucutre containing parameters used for clustering.
%
%
% By JMS, 4/14/2016
% Last edited 8/15/2017 for flexible mapping and masked sorting
%----------------------------------------------------------------------------------

% check inputs
p = check_inputs(varargin);

% get parameters
[npoints,nspikes] = size( spikes );
xvec = linspace( 0,npoints/Fs,npoints );

% set STATSET options for kmeans...use parallel computing if available
p.opt = statset( 'UseParallel',p.par );
    
% preallocate vector of valid spikes
kept = true( nspikes,1 );

%% EXTRACT WAVEFORM INFO
%==============================================

% first remove means of each spike
spikes = bsxfun( @minus,spikes,mean( spikes ) );

% Pull out features from the spikes dependent on the "method"
% implementation. For "raw", use components of the spike shape itself
if isnan( p.features ) && isnan( p.mapping )
    switch p.method
        case {'raw','RAW'}
            
            % Get the peak min/max of the spike
            [trough,peak] = spikeHeight( spikes,Fs );
            amp = peak - trough;

            % Get the AP slopes
            [dslope,uslope] = spikeHeight( diff( spikes ),Fs );

            % get the peak half-width
            halfwidth = halfWidth( spikes,trough,Fs );

            % area under curve (auc)
            auc = sum( spikes,1 );

            % energy
            energy = sum( spikes.^2 );

            % store into "features"
            features = zscore( [trough,peak,amp,dslope,uslope,...
                                halfwidth,auc',energy'] );
            Params.mapping = nan;

            % perform PCA if "decomp" is 1
            if p.decomp == 1
                [u,v] = svd( features' );
                p.level = min( p.level,size( features,2 ) );
                Params.mapping = u(:,1:p.level) * v(1:p.level,1:p.level);
                features = features * Params.mapping;
            end
            
        otherwise
            if strcmp( p.sortMethod,'mEM-GMM' ) || strcmp( p.sortMethod,'mVB-GMM' )
                [features,Params.mapping] = compute_spike_features( permute( spikes,[2,1,3] ),p.level,p.method,p.mask,p.distances,false );
            else
                [features,Params.mapping] = compute_spike_features( permute( spikes,[2,1,3] ),p.level,p.method,p.mask,p.distances,true );
            end
    end
elseif isnan( p.features ) && ~isnan( p.mapping )
    features = features * Param.mapping;
else
    features = p.features;
    Params.mapping = p.mapping;
end

% zscore the features for clustering stability
features = zscore( features ); 

%==============================================
%% Clustering
    
% perform the initial clustering
ID = perform_clustering();

% refine cluster using point-cluster probabilities
if all( isnan( p.mask ) )
    refine_cluster();  
end

% get the point-point distances for each cluster
if all( isnan( p.mask ) )
    get_density()
end

% store remaining parameters
Params.featureMethod = p.method;
Params.sortMethod = p.sortMethod;
Params.features = features;
Params.keptSpikes = kept;

% plotting
if p.plotting
    plotspikes();
end

%% Functions
    function p = check_inputs(inputs)
        % parse the optional inputs
        pnames = {'method','decomp','init','plotting','distances',...
            'par','search','level','reject','mapping','features','mask','sortMethod'};
        defaults = {'PCA',false,2,false,nan,false,false,5,.1,nan,nan,nan,'EM-GMM'};

        p = inputParser;             
        for j = 1:numel(pnames)
            p.addParameter( pnames{j},defaults{j} );
        end
        p.parse(inputs{:});
        p = p.Results;

        if (p.reject < 0) || (p.reject > 1)
            error( 'rejection probability not in the range [0,1]' );
        end
    end


    function ID = perform_clustering()
        % use a gaussian-mixture model EM algorithm to find the clusters
        
        % do our clustering using the p.init # of clusters or previous
        % clustering model. Store the new/updated model into Params
        switch p.sortMethod

            % === expectation maximization GMM ===
            case {'EM-GMM','em-gmm','EMGMM','emgmm','emGMM'} 
                if p.search && isnumeric( p.init )
                    maxclust = min( 6,floor( npoints*.1 ) );

                    % evaluate cluster number
                    E = evalclusters( features,'gmdistribution','Silhouette','klist',1:maxclust );
                    p.init = E.OptimalK;
                    fprintf( 'Optimal # of clusters: %i\n',p.init );
                end

                [~,Params.model,Params.likelihood,Params.prob]...
                    = mixGaussEm( features',p.init );

            % === variational bayes GMM ===
            case {'VB-GMM','vb-gmm','vb','VB'} 
                [~,Params.model,Params.likelihood] = mixGaussVb( features',p.init );
                Params.prob = Params.model.R;
                Params.model = rmfield( Params.model,'R' );
            
            % === masked EM-GMM ===
            case {'mEM-GMM','mem-gmm','maskedEM-GMM','maskedEM'}        
                if any( isnan( p.mask ) )
                    error( 'must provide mask matrix for masked EM-GMM' );
                end

                mask = zeros( size( features ) );
                for c = 1:size( p.mask,1 )
                    inds = c*p.level-p.level+1:c*p.level;
                    mask(:,inds) = repmat( p.mask(c,:)',1,p.level );
                end 

                [~,Params.model,Params.prob,Params.likelihood,Params.penalty]...
                    = maskedEM_GMM( features,mask,p.init,p.search );  

            % === masked VB-GMM ===
            case {'maskedVB','maskedvb','maksedVB-GMM','maskedvb-gmm','mVB-GMM','mvb-gmm'}
                % currently just a pseudo-masked VB algorithm using the "virtual ensemble" of X
                if any( isnan( p.mask ) )
                    error( 'must provide mask matrix for masked VB-GMM' );
                end
                
                mask = zeros( size( features ) );
                for c = 1:size( p.mask,1 )
                    inds = c*p.level-p.level+1:c*p.level;
                    mask(:,inds) = repmat( p.mask(c,:)',1,p.level );
                end
                
                X = compute_noise_ensemble( features,mask );
                [~,Params.model,Params.likelihood] = mixGaussVb( X.y',p.init );
                Params.prob = Params.model.R;
                Params.model = rmfield( Params.model,'R');

            % === k-means ===
            case {'k-Means','kmeans','KMEANS','K-Means','K-means'} 
                if p.search && isnumeric( p.init )
                    maxclust = min( 6,floor( npoints*.1 ) );

                    % evaluate cluster number
                    E = evalclusters( features,'kmeans','Silhouette','klist',1:maxclust );
                    p.init = E.OptimalK;
                    fprintf( 'Optimal # of clusters: %i\n',p.init );
                end

                Params.ID = kmeans( features,p.init );
        end

        % update the ID to remove small probabilities
        if isfield( Params,'prob' )
            Params.prob(:,sum(Params.prob) / size( features,1 )<0.01 ) = []; % remove any clusters containing only 1% of data
            [~,ID] = max( Params.prob,[],2 ); 
        end
    end


    function get_density()
        % finds the norm distance between pairs of points in each cluster, then
        % sums these distances and divides by N^2, where N = number of points
        % in cluster c. Finally sums the normalized values together into W
        
        uID = unique( ID );
        uID(uID==0) = []; % remove IDs = 0
        nclust = numel( uID );
        Params.clusterDensity = zeros( nclust,nclust );
        
        % compute the pairwire distances
        D = pdist2( features,features );
               
        % loop over the individual clusters and store the average density
        % between points from clusters i and j
        for i = uID
            for j = uID
                distance = range( mean( D(ID==i,ID==j),2 ) );
                [~,volume] = convhulln( features(ID==i | ID==j,:) );
                Params.clusterDensity(i,j) = distance / volume;
            end
        end
    end


    function refine_cluster()
        % refine each cluster to discard points with less than a specific probability
        % of belonging to that cluster. This increases the specificity of
        % detected spikes at the cost of the total number of kept spikes...this
        % can be troublesome if the clusters are largely overlapping, as this
        % will throw away most of the points that exist within the overlap.

        % parse the probability matrix into each cluster
        nClust = numel( unique( ID ) );
        for i = 1:nClust
            if nClust == 1 % if only one group specified
                bad = medoutlier( Params.prob,4 ); % remove outliers
                ID(bad) = 0;
                kept(bad) = 0;
            else
                index = ID==i;
                if sum( index ) <= p.level
                    ID(index) = 0; % make these equal to zero;
                else
                    % if probability of being in cluster "i" is less than
                    % p.reject, drop the point
                    bad = Params.prob(index,i) <= p.reject; 
                    ID(index(bad)) = 0; % change ID to 0 for these bad spikes
                    kept(index(bad)) = 0; % change kept for these spikes to 0
                    clear index bad 
                end
            end
        end
    end        


    function plotspikes()
        % plot the spike waveforms and feature space based on color.
        % Gray waveforms/feature points represents dropped spikes (ID == 0)
        % Additionally, the size of the feature points represents the
        % probability of belonging to the nearest centroid

        uID = unique( ID );      
        uID(isnan( uID )) = [];
        figure; 
        col = {[.4 .6 1],[.5 0 .6],[.9 .6 0],[.2 .7 .2],[.9 .2 0],[0 .6 .7]};

        % get the Area matrix for the scatter plot
        Area = Params.prob.*100;
        Area = Area.^2;
        Area = bsxfun( @rdivide,Area,max( Area ) ) .* 5;

        % first plot waveforms/features in gray for ID == 0
        xvec = linspace( 0,size(spikes,1)/Fs,size(spikes,1) );
        if sum(uID==0) >= 1
            subplot(2,1,1);
                plot(xvec,spikes(:,ID==0),'color',[.8 .8 .8]);

            subplot(2,1,2);
                s = scatter3( Params.features(ID==0,1),Params.features(ID==0,2),...
                    Params.features(ID==0,3),max( Area(ID==0) ).^2' );
                xlabel('feature 1'); ylabel('feature 2'); zlabel('feature 3');
                set(s,'markerfacecolor',[.8 .8 .8],'markeredgecolor','none');
            nID = numel(uID)-1;
        else
            nID = numel(uID);
        end

        % now plot the sorted spikes for each ID > 0
        for j = 1:nID
            subplot(2,1,1); hold on;
                plot( xvec,spikes(:,ID==j),'color',col{j} );

            subplot(2,1,2); hold on;
                s = scatter3( Params.features(ID==j,1),Params.features(ID==j,2),...
                    Params.features(ID==j,3),Area(ID==j,j).^2 );
                xlabel('feature 1'); ylabel('feature 2'); zlabel('feature 3');
                set(s,'markerfacecolor',col{j},'markeredgecolor','none');
        end
        
        % clean up the tick labels
        subplot(2,1,1); set(gca,'tickdir','out','box','off');
        subplot(2,1,2); set(gca,'tickdir','out','box','off');      
    end
end