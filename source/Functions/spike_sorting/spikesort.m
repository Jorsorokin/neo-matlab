function [ID,Params] = spikesort(spikes,Fs,varargin)
% ------------ [ID,Params] = spikesort(spikes,Fs,varargin) ---------
%
%   Sort columns extracellular spikes using K-means clustering by pulling
%   out relevant features including peak-to-trough height, spike
%   half-width, AHP, etc. Optionally uses kernel PCA with a gaussian kernel
%   for dimensionality reduction before clustering. In addition the function 
%   can iteratively sort using 1:20 clusters and approximate the true number
%   of clusters in the dataset.
%
%   This function can also work in a parallel computing environment.
%
%                   >>> INPUTS >>>
%
% Required:
%   spikes: m x n matrix of spike snips with m = length of each spike waveform,
%               n = # of spikes
%   Fs: the sampling rate (used for upsampling and width calculation)
%
% Optional:
%           * All optional inputs should be given as name-value pairs *
%
%   method: string specifying the feature-selection method to use for
%           sorting. Valid arguments are { "pca", "ica", "raw" }. For
%           "raw", measurements from the actual waveforms such as spike
%           height, width, etc. will be used. 
%
%   decomp: boolean (0 or 1)...if 1, will perform dimensionality reduction
%           via kPCA on the extract waveform features (amplitude, ahp, etc). 
%           Note, this only applies if "raw" is chosen as "method"  
%   
%   init: the number of clusters for EM algorithm or a model structure from
%         a previous run of an EM algorithm (containing the mu, sigma, and 
%         weights for each cluster). Default is to set to 2 clusters and
%         not use a previous model.
%   
%   PCs: an n x m matrix of principal components to be projected onto.
%        Only applies if "method" equals "pca"
%
%   features: an n x m matrix of features computed from a previous spike
%             sorting run. If provided, no PCs will be computed
%
%   usamp: boolean 0 or 1...if 1, will upsample by x4 (default = 0)
%
%   level: the total number of featuresto keep for clustering. By default,
%          an automated feature-selection algorithm will determine the
%          relevant features. If "level" is specified, the algorithm will
%          keep the top "level" number of features.
%
%   reject: the min probability of point X belonging to any cluster.
%           probabilities lower than this will result in rejection of point
%           X from the cluster analysis (i.e. ID = 0). Default value is 0.6
%           (can be any value in [0,.01,.02 ... .99, 1.0]) 
%   
%   plotting: 0 or 1...plots the sorted spike snips and scatter plot of
%             first three features or PC space (default = 0).
%             Additionally plots the within/between distance and R value if
%             "search" is set to 1. 
%   
%   par: 0 or 1...if 1, sets the "UseParallel" option to true for Kmeans
%        and invokes parallel computing environment
%   
%   search: 0 or 1...if 1, searches for an optimal cluster number using a
%           GMM distribution
%
%                   <<< OUTPUTS <<<
%
% ID: n x 1 vector of cluster ID's 
%
% Params: strucutre containing parameters used for clustering.
%
%
% By JMS, 4/14/2016
%----------------------------------------------------------------------------------

% check inputs
p = check_inputs(varargin);

% get parameters
[npoints,nspikes] = size( spikes );
xvec = linspace( 0,npoints/Fs,npoints );
if p.usamp == 1
    xxvec = resample( xvec,4,1 ); % 4x up-sample rate
    Fs = Fs*4; 
    spikes = spline(xvec,spikes',xxvec)'; 
    xvec = xxvec;
end

% set STATSET options for kmeans...use parallel computing if available
p.opt = statset('UseParallel',p.par);
    
% preallocate vector of valid spikes
kept = true( nspikes,1 );

%% EXTRACT WAVEFORM INFO
%==============================================

% first remove means of each spike
spikes = bsxfun( @minus,spikes,mean( spikes ) );

% Pull out features from the spikes dependent on the "method"
% implementation. For "raw", use components of the spike shape itself
switch p.method
    case {'raw'}
        
        if isnan( p.features )
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

            % perform PCA if "decomp" is 1
            if p.decomp == 1
                if p.PCs == 0
                    [u,v] = svd( features' );
                    p.level = min( p.level,size( features,2 ) );
                    Params.PCs = u(:,1:p.level) * v(1:p.level,1:p.level);
                else
                    Params.PCs = p.PCs;
                end
                features = features * Params.PCs;
            end
        end
            
    case {'pca'}
        
        % check if "p.PCs" is supplied. If so, project onto the provided
        % PCs, else perform PCA
        
        if isnan( p.features )
            if p.PCs == 0
                % perform PCA on the spike waveforms
                p.level = min( p.level,nspikes );
                [u,v] = svd( spikes );
                eigval = diag( v ).^2;

                % only keep those that explain 95% of the variance
                eigsum = eigval / sum( eigval );
                eigsum = cumsum( eigsum );
                ind = find( eigsum>=.95,1 );
                if ~isempty( ind ) && ind < p.level
                    p.level = ind;
                end

                % create our PCs
                Params.mapping = u(:,1:p.level) * v(1:p.level,1:p.level);
                Params.eigval = eigval;
            else
                Params.mapping = p.PCs;
            end
        
            % project onto the PCs
            if isnan( p.features )
                features = spikes' * Params.mapping;
                p.level = size( Params.mapping,2 );
            end
        end

    case {'ica'}
        
        if isnan( p.features )       
            % perform ICA on the spike waveforms
            [Params.mapping,features,~] = fastica( spikes','numOfIC',p.level,...
                'lasteig',min( p.level*2,size(spikes,2) ),...
                'Verbose','off','stabilization','on' );
        end
    
    otherwise
        if isnan( p.features )
            [features,Params.mapping] = compute_mapping( spikes',p.method,p.level ); % use defaults
        end         
end

% zscore the features for clustering stability
if isnan( p.features )
    features = zscore( features ); 
else
    features = p.features;
end
%==============================================

%% Clustering
    
% perform the initial clustering
ID = perform_clustering();

% refine cluster using point-cluster probabilities
refine_cluster();  

% get the point-point distances for each cluster
get_density()

% store remaining parameters
Params.featureMethod = p.method;
Params.features = features;
Params.keptSpikes = kept;

% plotting
if p.plotting == 1
    plotspikes();
end

%% Functions
    function p = check_inputs(inputs)
        % parse the optional inputs
        pnames = {'method','decomp','init','usamp','plotting',...
            'par','search','level','reject','PCs','features'};
        defaults = {'raw',0,2,0,0,0,0,5,.6,0,0};
        options = {{'raw','ica','pca'},{0,1},{nan},{0,1},{0,1},...
                         {0,1},{0,1},{1:50},{linspace(0,1,101)},{nan},{nan}};

        p = inputParser;             
        % loop over the rest of the optional inputs
        for j = 1:numel(pnames)
            if ischar(options{j}{1})
                p.addParameter( pnames{j},defaults{j},@(x) max(strcmp(x,options{j})) == 1 );
            else
                if strcmp( pnames{j},'init' ) || strcmp( pnames{j},'PCs' ) || strcmp( pnames{j},'features' )
                    p.addParameter( pnames{j},defaults{j},@(x) ~isempty(x) );
                else
                    p.addParameter( pnames{j},defaults{j},@(x) max(x == cell2mat(options{j})) == 1 );
                end
            end
        end

        p.parse(inputs{:});
        p = p.Results;
    end

    function ID = perform_clustering()
        % use a gaussian-mixture model EM algorithm to find the clusters
        
        % only perform search if previous model not supplied
        if p.search == 1 && isnumeric( p.init )        
            maxclust = min( 6,floor( npoints*.1 ) );
            
            % evaluate cluster number
            E = evalclusters( features,'gmdistribution','Silhouette','klist',1:maxclust );
            p.init = E.OptimalK;
            fprintf( 'Optimal # of clusters: %i\n',p.init );
        end
        
        % do our clustering using the p.init # of clusters or previous
        % clustering model. Store the new/updated model into Params
        [ID,Params.model,~,Params.prob] = mixGaussEm( features',p.init );
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