function [labels,model,prob,likelihood] = tDistMixEM( X )
    % [labels,model,prob,likelihood] = tDistMixEM( X )
    %
    % Clusters points in X using the Expectation Maximization for T-distributed Mixture Modeling
    %
    % Inputs:
    %   X - an n x d matrix of points to be clustered
    %
    % Outputs:
    %   labels - an n x 1 vector of cluster assignment for each point 
    %
    %   model - a structure containing the fitted params:
    %       mu - cluster means
    %       sigma - cluster covariances
    %       dof - the DOF of the t-dist
    %       weights - cluster weights
    %
    %   prob - an n x K matrix of probabilities for each point 1:n and cluster 1:K
    %
    %   likelihood - the log likelihood function for each iteration. The final clustering 
    %                is associated with the maximum of the likelihood 
    %
    %
    % References:
    %   Shohan et al. 2002, J. Neuro. Methods
    %
    % Written by Jordan Sorokin
    % 10/7/2017

    %% Initialization
    [n,m] = size( X );
    Kmin = 1;                           % min cluster numbers
    K = 100;                            % default to 100 clusters to begin 
    labels = kmeans( X,Kinit );         % first pass clustering 
    p = zeros( 1,K ) + (1 / K);         % evently distributed priors
    mu = get_clust_params( X,labels );  % initial means
    sigma = eye( m,m,K );               % initial covariances = I
    dof = 50;                           % initial DOF 
    Lmax = -inf;                        % initial likelihood
    minErr = 1e-6;                      % for convergence

    % create the model. We will update the model fields when 
    % the likelihood L exceeds the current likelihood Lmax
    model = struct( 'mu',mu,'sigma',sigma,'weights',p,'dof',dof );

    %% Main Run
    while (K >= Kmin) && (err >= minErr)

        % E-STEP


        % M-STEP


        % CHECK CONVERGENCE
        if L > Lmax
            model.mu = mu;
            model.sigma = sigma;
            model.weights = p;
            model.dof = dof;

            Lmax = L; % new likelihood maximum
            K = K - 1; % update K 
        end
    end



    function [mu,sigma] = get_clust_params( X,labels )
        uID = unique( labels )';
        nK = numel( uID );
        mu = zeros( nK,: );
        sigma = zeros( m,m,nK );

        % compute means & covariances
        for j = uID
            thesePts = X(labels==j,:);
            npt = size( thesePts,1 );
            mu(j,:) = mean( thesePts );
            sigma(:,:,j) = (thesePts' * thesePts) / npt;
        end
    end



end