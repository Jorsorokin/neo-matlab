function [label, model, R, likelihood, penalty, splitClusters, killedClusters] = maskedEM_GMM( X,mask,k,varargin )
    % [label, model, likelihood, R, penalty,splitClusters, killedClusters]
    %       = maskedEM_GMM( X,mask,k,maxSplitAttempts ) 
    %
    % Perform EM algorithm on a masked data matrix X for fitting a 
    % gaussian mixture model with a predefined number of components
    %
    %                   >>> INPUTS >>> 
    %   X: 
    %       n x d feature matrix, with n = observations, d = variables
    %
    %   mask: 
    %       n x d sparse mask matrix, with each element (i,j) in the set [0,1]
    %
    %   k:
    %       the number of clusters 
    %
    %   (maxSplitAttempts):
    %       scalar determining how many attempts per cluster the algorithm should 
    %       try to split that cluster to improve the likelihood function. Default
    %       = 0 (i.e. no splitting attempts)
    %
    %                   <<< OUTPUTS <<<
    %   label: 
    %       n x 1 cluster label for each ith data point
    %
    %   model: 
    %       trained model structure with fields:
    %           mu - d x k matrix of means of each cluster
    %           sigma - d x d x k tensor of covariances of each cluster
    %           w - 1 x k vector weights of each cluster
    %
    %   R: 
    %       n x k matrix of probabilities of spike i belonging to cluster c
    %
    %   likelihood: 
    %       vector of the full log likelihood for each iteration j
    %
    %   penalty:
    %       vector of the modified AIC penalty function for each iteration j
    %
    %   splitClusters:
    %       boolean vector indicating if any clusters were split on each iteration j
    %
    %   killedClusters:
    %       boolean vector indicating if any clusters were killed on each iteration j
    %   
    % 
    % original mixGaussEm by Mo Chen
    %
    % adapted for masked EM by JMS, 8/8/2017

    % globals
    tol = 1e-8;
    iter = 1;
    maxiter = 500;
    likelihood = -inf( 1,maxiter );
    penalty = -inf( 1,maxiter );
    killedClusters = zeros( 1,maxiter );
    splitClusters = zeros( 1,maxiter ); 
    
    % splitting parameters
    splitAttempts = zeros( 1,k );
    if nargin > 3 && ~isempty( varargin{1} )
        maxSplitAttempts = varargin{1};
    else
        maxSplitAttempts = 0;
    end

    % error checking
    [n,d] = size( X );
    if ~all( size( mask ) == [n,d] )
        error( 'Mask and data matrices not equal sizes' );
    end
    if k > d
        warning( ['Number of requested clusters is greater than dimension of X',...
            'Reducing number of clusters to equal d-1'] );
        k = d-1;
    end

    % compute our noise ensemble given the data and mask
    virtualX = compute_noise_ensemble( X,mask );
            
    % and +1 component for noise and initialize random labels 
    k = k + 1; 
    noise = randn( n,d );
    noise_mu = mean( noise )';
    noise_sigma = cov( noise );
    [label,R] = initialization( k,n ); % our first "expectation"
    clear noise

    % MAIN LOOP
    % ================================================================
    while iter <= maxiter 
        iter = iter + 1;
        numDeleted = 0; % global
        numSplit = 0; % global

        % run one iteration of the EM algorithm for k gaussian components
        % and automatically kill any clusters with singular covariances
        [label,R,badClust] = reassign_points( label,R ); 
        [label,model,R,lh] = EM_one_step( R,label );
        killedClusters(iter) = numDeleted;
        
        % compute the penalty function
        penalty(iter) = penalization( label,lh );

        % compute change in likelihood
        %likelihood(iter) = sum( lh ); % full likelihood function
        dLh = abs( lh - likelihood(iter-1) ) / abs( lh );
        
        % check if we should try splitting clusters
        splitAttempts( badClust ) = [];
        if (dLh <= tol) && (maxSplitAttempts > 0) && any( splitAttempts ~= maxSplitAttempts )
            try_splitting_clusters(); % recursively calls EM_GMM() for each cluster separately, changes R, model, and lh in place
            splitClusters(iter) = numSplit;
            dLh = abs( lh - likelihood(iter-1) ) / abs( lh );
            penalty(iter) = penalization( label,lh );
        end
        
        % store the likelihood
        likelihood(iter) = lh;
        if penalty(iter) > max( penalty(1:iter-1) )
            bestModel = model;
            bestR = R;
            bestIter = iter;
        end

        % now check for convergence
        if dLh <= tol
            break
        end
    end
    % ================================================================
    
    if iter > maxiter
        warning( 'Reached maximum # of iterations before converging' );
        iter = maxiter;
    end

    % take the best model possible as our solution
    %model = bestModel;
    %R = bestR;
    [~,label] = max( R,[],2 );
    bestIter = iter;
    likelihood = likelihood(2:bestIter);
    penalty = penalty(2:bestIter);
    splitClusters = splitClusters(2:bestIter);
    killedClusters = killedClusters(2:bestIter);


    %% FUNCTIONS
    function [label,model,R,lh] = EM_one_step( R,label )
        % runs one step of the EM algorithm given the current probabilities and labels
        % and outputs new probabilities, labels, mixture model, and likelihood    
        
        % ----- M STEP -----
        % compute the model given current cluster assignments 
        model = maximization( R,label );
        % ------------------

        % ----- E STEP ----- 
        % compute the posterior probabilities given the cluster model
        R = expectation( R,virtualX.y,virtualX.eta,model );
        % ------------------
        
        % compute the new likelihood and labels
        [R,lh] = compute_likelihood( R,model.w );
        [~,label(:)] = max( R,[],2 );
    end

    function [label,R] = initialization( init,n )
        % initializes the label vector and probability matrix 
        
        if isscalar( init ) % number of clusters
            label = ceil( init*rand( n,1 ) );  
        elseif isvector( init ) && (numel( init ) == n) % cluster assignment
            label = init;
        end
        R = full( sparse( 1:n,label,1 ) ); 
    end

    function R = expectation( R,X,eta,model )
        % computes the log gaussian pdf for each point and each cluster 
        
        [~,nClust] = size( R );
        for i = 1:nClust
            
            [U,singular] = chol( model.Sigma(:,:,i) ); % upper triangular matrix
            if singular ~= 0
                R(:,i) = 0;
                continue
            end

            Q = bsxfun( @minus,X,model.mu(:,i)' )/U;
            invU = inv( U ); 
            
            R(:,i) = -0.5 * (...
                d*log( 2*pi )...
                + 2*log( sum( diag( U ) ) )... % equivalent to log( det( Sigma ) );
                + dot( Q,Q,2 )... % equivalent to: X * inv( Sigma) * X'
                + sum( bsxfun( @times,eta,diag( invU * invU' )' ),2 )...
                ); 
        end
    end

    function model = maximization( R,label )    
        % Computes the current model of the clusters 1:k given the data points for each
        
        % reassign any points in which only one point is assigned to a cluster
        [nPoints,nClust] = size( R );
        [clustSize,clustIDX] = count_cluster_points( label );

        mu = zeros( d,nClust );
        Sigma = zeros( d,d,nClust );
        mask_intersect = zeros( d,d );

        % now loop over non-empty clusters and compute covariances
        for i = 1:nClust-1
            nz = clustSize(i);
            C = clustIDX(:,i);
            L = mask(C,:) > 0; % all unmasked points
            M = ~L; % all masked points
            doubleM = M;

            mu(:,i) = (sum( virtualX.y(C,:).*L ) + sum( M ).*virtualX.mu) / nz; % the cluster means
            eta = sum( virtualX.eta(C,:).*L ) + sum( M ).*virtualX.sigma; % variances of masked features
            yhat = virtualX.y(C,:) - mu(:,i)'; % de-meaned virtual ensemble
            
            for j = 1:d
                %yhat = virtualX.y(C,:) - mu(j,i); % subtracts the jth feature mean from all features in y
                doubleM(:) = (M(:,j) & M); % where M_i intersect M_j (both masked)
                mask_intersect(j,:) = sum( doubleM ); % |M_i intersect M_j|
                Sigma(j,:,i) = yhat(:,j)' * (yhat .* ~doubleM);
            end

            muhat = virtualX.mu' - mu(:,i);
            Sigma(:,:,i) = (Sigma(:,:,i)...% covariance of unmasked features
                            + mask_intersect .* (muhat * muhat')... % covariance of masked features
                            + diag( eta )) / nz; 
        end
        Sigma(:,:,end) = noise_sigma;
        mu(:,end) = noise_mu;
        model = struct;
        model.mu = mu;
        model.Sigma = Sigma;
        model.w = sum( R ) / nPoints;
        model.w(end) = min( model.w(end),0.05 ); % small constant prior ensured
    end

    function penalty = penalization( label,lh )
        % calculates the model penalty using a modified AIC
        
        freeParams = count_free_parameters( label );
        penalty = 2 * (sum( freeParams(1:end-1) )-1 - sum( lh ) ); % AIC
    end

    function [R,lh] = compute_likelihood( R,weights )
        % computes the likelihood for each cluster k as:
        %   SUM_1:N{ log{ P(Y=k | x,theta) * P(Y == k) } } / N
        
        nPoints = size( R,1 );
        R = bsxfun( @plus,R,log( weights ) ); % log( P(Y=k | x,theta) ) + log( P( theta ) ) = log{ P(Y==k | x,theta) * P(theta) }
        T = logsumexp( R,2 );
        R = exp( bsxfun( @minus,R,T ) ); % log( P(Y==k |x,theta)*log( P(theta) / SUM_K{ P(Y==k |x,theta)*log( P(theta) } )  )
        lh = sum( T ) / nPoints; 
    end

    function freeParams = count_free_parameters( label )
        % calculates the number of free parameters for each component

        [clustSize,clustIDX] = count_cluster_points( label ); 
        nClust = numel( clustSize );
        freeParams = zeros( 1,nClust );
        for i = 1:nClust
            r = sum( mask(clustIDX(:,i),:),2 ); % sum of unmasked features for each data point in this cluster
            freeParams(i) = sum( (r .* (r+1)) / 2 + r + 1 ) / clustSize(i);
        end
    end  

    function [label,R,badClusters] = reassign_points( label,R )
        % reassign points away from any cluster that has less than 3 points

        [clustSize,clustIDX] = count_cluster_points( label );
        badClusters = find( (clustSize(1:end-1) <= 1 ) );
        if ~isempty( badClusters )
            for j = badClusters
                badPts = clustIDX(:,j);
                R(badPts,:) = 0;
                R(badPts,end) = 1;
            end
            R(:,badClusters) = [];
            [~,label(:)] = max( R,[],2 );
            numDeleted = numDeleted + numel( badClusters );
            k = k - numDeleted;
        end
    end

    function [clusterSize,clusterIDX] = count_cluster_points( label )
        % returns the number of points in each cluster and the indicies of the points in each cluster

        clusterIDX = bsxfun( @minus,repmat( 1:k,n,1 ),label ) == 0; % finds which indices of "label" belong to each cluster
        clusterSize = sum( clusterIDX );
    end

    function try_splitting_clusters()
        % recursively call maskedEM_GMM with subsets of data for each identified cluster, and split each cluster
        % if it improves the likelihood of that cluster or reduces the bi/multi-modal distribution of the cluster

        [clustSize,clustIDX] = count_cluster_points( label );

        fprintf( 'Attempting to increase the log likelihood by splitting clusters' );
        %oldModel = struct;
        for i = 1:k-1
            fprintf( '.' );

            % update the splitting history for this cluster, even if the subsequent split attempt fails. 
            % This prevents the algorithm from attempting to split any ith cluster infinitely many times
            splitAttempts(i) = splitAttempts(i) + 1;

            % check if we should attempt to split this cluster
            if splitAttempts(i) == maxSplitAttempts || clustSize(i) <= 8 % splitting will result in subsequent killing of at least one cluster, so just ignore
                continue 
            end

            % compute the old likelihood of this cluster, only considering the points assigned to it
            %   SUM{ log{ P( X_k == k | theta_k ) * P( theta_k ) } }, X_k denotes sub-sample 
               
            % oldModel.mu = model.mu(:,i);
            % oldModel.Sigma = model.Sigma(:,:,i);
            % oldModel.w = model.w(i);
            % oldProb = zeros( clustSize(i),1 );
            % oldProb = expectation( oldProb,virtualX.y(C,:),virtualX.eta(C,:),oldModel );
            % [~,oldLh] = compute_likelihood( oldProb,model.w(i) ); 
            
            % call maskedEM_GMM() with the subset of data and no_splitting
            C = clustIDX(:,i); 
            [~,model2,~,~] = maskedEM_GMM( X(C,:),mask(C,:),2,0 ); % 2 clusters
            if numel( model2.w ) < 2
                continue % the two clusters were re-merged automatically, so don't do anything else  with this cluster
            end
            
            % normalizes the weights of the splitted clusters and data point probabilities: 
            %   SUM{ P(theta_s) } = P(theta_k), for s = [1,2]
            newModel.w = model2.w * model.w(i); 
            
            % remove noise terms
            newModel = model;
            newModel.w(end) = [];
            newModel.Sigma(:,:,end) = [];
            newModel.mu(:,end) = []; 
            
            % add the model from the sub EM run
            newModel.mu(:,i) = model2.mu(:,1);
            newModel.mu(:,end+1) = model2.mu(:,2);
            newModel.Sigma(:,:,i) = model2.Sigma(:,:,1);
            newModel.Sigma(:,:,end+1) = model2.Sigma(:,:,2);
            newModel.w(i) = model2.w(1);
            newModel.w(end+1) = model2.w(2);
            
            % re-add the noise terms
            newModel.w(end+1) = model.w(end);
            newModel.Sigma(:,:,end+1) = model.Sigma(:,:,end);
            newModel.mu(:,end+1) = model.mu(:,end);
            
            % compute the probabilities
            newProb = zeros( n,numel( newModel.w ) );
            newProb = expectation( newProb,virtualX.y,virtualX.eta,newModel ); % P( X_k = s | theta_s ) normalized by P(theta_k) from above
            newLabel = max( newProb,[],2 );
            newModel = maximization( newProb,newLabel );
            if any( all( newProb==0 ) )
                continue % at least one cluster Sigma was rank deficient (singular assignment)
            end

            % compute new likelihood, and compare with previous likelihood
            %   SUM_s{ SUM_N{ log{ P( X_k = s | theta_s ) * P( theta_s ) } } } 
            %   vs. SUM_N{ log{ P( X_k = k | theta_k ) * P( theta_k } }
            [newProb,newLh] = compute_likelihood( newProb,newModel.w );
            if newLh > lh
                k = k + 1;
                numSplit = numSplit + 1;                
                model = newModel;

                % alter the current cluster model and add a new one
                % model.mu(:,i) = newModel.mu(:,1); 
                % model.mu(:,k) = newModel.mu(:,2);
                % model.Sigma(:,:,i) = newModel.Sigma(:,:,1);
                % model.Sigma(:,:,k) = newModel.Sigma(:,:,2);
                % model.w(i) = newModel.w(1);
                % model.w(k) = newModel.w(2);

                % updates the child cluster's split histories to match the parent
                splitAttempts(k) = splitAttempts(i);
            end
        end

        if numSplit > 0
            fprintf( '\nSuccessfully split %i clusters\n',numSplit );
            %R = zeros( n,k );
            %R = expectation( R,virtualX.y,virtualX.eta,model );
            %[R,lh] = compute_likelihood( R,model.w );
            R = newProb;
            lh = newLh;
            [~,label] = max( R,[],2 );
        end
    end
end 