function [L,W,Y,model,P,dRatio] = dslClust( X,varargin )
%
% [L,W,Y,model,P,dRatio] = dslClust( X,(W0,k0,ndim,search) )
% 
% clustering via Discriminant Subspace Learning (DSL)
%
% uses linear discriminant analysis (LDA) and the expectation maximization 
% (EM) algorithm to cluster the data in the matrix X into distinct groups, 
% assuming a Gaussian Mixture Model (GMM) prior. 
%       * For details, see Kesthkaran & Yang 2017
%
% Inputs:
%   X - an n x m matrix, with n = number of samples per observation, m =
%       number of observations.
%   (W) - the initilization projection matrix. If not provided, W0 will
%          be calculated via LDA-kmeans
%   (k0) - the initial cluster number. Default = 1
%   (ndim) - the number of dimensions in the reduced subspace (default = 3)
%   (search) - boolean flag for searching for the optimal cluster number.
%              If 0, the function will perform one iteration of LDA-GMM
%              using k0 # of clusters; else, it will sweep through a
%              variety of cluster numbers and automatically detect the most
%              appropriate (see Kesthkaran & Yang paper). Default = true.
%
% Outputs:
%   L - a m-component vector with cluster assignments (1 -> k) for each
%       observation in X
%   W - the final projection matrix used for projecting onto the best subspace
%   Y - the projections of the spikes onto the final W projection matrix
%   model - a structure containing the means (mu), covariances (sigma), and
%           probability of each cluster (w)
%   P - the probabilities of each point belonging to each cluster     
%   dRatio - the ratio of within- and between-distances of the cluster
%            centroids
%
% By JMS, 7/17/2017

% initialize the parameters
MAXITER = 100;
MINTOL = 1e-4;
[n,m] = size( X );

% check inputs
if nargin > 1 && ~isempty( varargin{1} )
    W = varargin{1};
else
    W = [];
end
if nargin > 2 && ~isempty( varargin{2} )
    k0 = varargin{2};
else
    k0 = 1;
end
if nargin > 3 && ~isempty( varargin{3} )
    ndim = varargin{3};
else
    ndim = 3;
end
if nargin > 4 && ~isempty( varargin{4} )
    search = varargin{4};
else
    search = true;
end

% find W0 if necessary
if isempty( W )
    %W = initialize_W();
    [u,s] = svd( X );
    W = u(:,1:ndim) * s(1:ndim,1:ndim);
end

% Loop over the cluster numbers for automatic cluster identification
if search
    K = k0 + 1;
    Khat = K;
    fprintf( 'Searching for optimal cluster #' );
    while K <= 15 % arbitrarily chose 15 here
        fprintf( ' . ' );

        % run LDA-GMM assuming k2 number of clusters
        [~,W] = lda_gmm( W,X,K );

        % now run k-means to cluster the data
        Y = W' * X;
        % whiten Y
        [labels,mu] = kmeans( Y',K );

        % loop over the cluster numbers, compute the distribution of
        % projections of clusters onto cluster-to-cluster vectors. 
        for i = 1:K    
            % get samples X_i (those in cluster i)
            Y_i = Y(:,labels==i);

            for j = 1:i-1   
                % get samples X_j
                Y_j = Y(:,labels==j);

                % project onto vector mu_i - mu_j
                muVec = mu(i,:) - mu(j,:);
                Yproj_i = muVec * Y_i;
                Yproj_j = muVec * Y_j;

                % compute the histograms of each X_i, X_j, and the two combined
                % and the anderson-darling test statistic
                if i ~= j
                    [~,~,adstat_i] = adtest( Yproj_i );
                    [~,~,adstat_j] = adtest( Yproj_j );
                    [~,~,adstat_ij] = adtest( [Yproj_i,Yproj_j] );

                    % compare the statistics. If any individual X_i or X_j has a
                    % smaller statistic than the combined X_ij, then we have too
                    % many clusters
                    if (adstat_i > adstat_ij) && (adstat_j > adstat_ij)
                        Khat = Khat - 1;
                    end    
                end
            end
        end
        % now check if Khat < K. If so, we have surpased the number of
        % clusters that actually exist in our data
        if Khat < K
            break;
        else
            K = K+1;
            Khat = K;
        end
    end
    k0 = Khat;
    fprintf( '\nOptimal # of clusters: %i\n',k0 );
else
    fprintf( 'Clustering via LDA-GMM...\n' );
    [~,W] = lda_gmm( W,X,k0 );
end

% finally, compute our final clustering 
Y = W' * X;
Y = whiten_data( Y );
[L,model,~,P] = mixGaussEm( Y,k0 );
dRatio = get_dratio( X,L,W );

    %% Helper functions
    function W = initialize_W()
        % initialize the projection matrix W   
        i = 1;
        msub = ceil( 0.1*m ); % 10%
        fprintf( 'Initializing projection matrix W...\n' );
        while (i <= 10) % arbitrarily set to 10 iterations
            
            % extract a random subset of X
            subset = randperm( m,msub ); % 10% of observations
            X_sub = X(:,subset);

            % perform PCA to get the projection matrix W = eigenvectors of X
            [u,s] = svd( X_sub );
            W0 = u(:,1:ndim) * s(1:ndim,1:ndim);

            % perform LDA-kmeans to find true W0
            [labels,W0,mu] = lda_km( W0,X_sub,k0 );
            
            % find the within/between-distance ratio
            if i > 1
                [dr,idx] = max( [dr,get_dratio( X_sub,labels,W0 )] );
                if idx == 2
                    W = W0;
                end
            else
                dr = get_dratio( X_sub,labels,W0 );
                W = W0; % make a copy;
            end
            i = i+1;
        end
    end

    function [labels,W,mu] = lda_km( W,data,k )
        % use linear discriminant analysis (LDA) and k-means to cluster
        % data into distinct clusters
        i2 = 1;
        while (i2 <= MAXITER)
            
            % compute projections
            Y = W' * data;
            
            % whiten the projections
            Y = whiten_data( Y );
            
            % run k-means
            [newlabels,mu] = kmeans( Y',k );
            
            % run LDA
            [~,W] = gda( data',newlabels,k-1 );
            W = W.M;
            
            % get the distance ratio
            newRatio = get_dratio( data,newlabels,W );
        
            % compare with previous labels
            if i2 > 1
                if abs( oldRatio - newRatio ) <= MINTOL
                    break;
                end
            else
                labels = newlabels;
                oldRatio = newRatio;  
            end
            i2 = i2+1;
        end 
    end

    function [labels,W] = lda_gmm( W,data,k )
        % use linear discriminant analysis (LDA) and expectation maximization
        % to cluster data into distinct clusters
        i3 = 1;
        while (i3 <= MAXITER)
            
            % compute projections
            Y = W' * data;
            
            % whiten the projections
            Y = whiten_data( Y );
            
            % run k-means
            [newlabels,model] = mixGaussEm( Y,k );
            
            % run LDA
            [~,W] = lda( data',newlabels,k-1 );
            W = W.M;
            
            % get the new distance ratio
            newRatio = get_dratio( data,newlabels,W );
        
            % compare with previous ratio
            if i3 > 1
                if abs( oldRatio - newRatio ) <= MINTOL
                    break;
                end
            else
                labels = newlabels;
                oldRatio = newRatio;       
            end
            i3 = i3+1;
        end 
    end

    function dRatio = get_dratio( data,labels,W )
        % find the ratio of between-/within-cluster distances
        
        % create our label matrix
        [ndims,npts] = size( data );
        nclust = max( labels );
        labelMat = zeros( npts,nclust );
        M = zeros( ndims,nclust );
        for k = 1:nclust
            labelMat(labels==k,k) = 1;
            M(:,k) = mean( data(:,labels==k),2 );
        end
        
        % compute the within/between distances
        maskedMu = M * labelMat';
        withinDist = (data - maskedMu) * (data - maskedMu)';
        betweenDist = maskedMu * labelMat * M';
        
        % compute the ratio
        dRatio = trace( (W' * betweenDist * W) / (W' * withinDist * W) );
    end
    
    function Y = whiten_data( Y )
        % whiten the input data
        C = cov( Y' );
        U = C^(-1/2);
        Y = U * Y;
    end

end