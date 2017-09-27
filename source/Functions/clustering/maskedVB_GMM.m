function [model, likelihood] = maskedVB_GMM(X, m, mask, prior )
    % Variational Bayesian inference for Gaussian mixture.
    %
    % Inputs: 
    %   X: d x n data matrix
    %   m: k (1 x 1) or label (1 x n, 1 <= label(i) <= k) or model structure
    %   mask: a d x n masking matrix (0 <= (i,j) <= 1)
    %   (prior): structure of priors over hyperparameters (optional)
    %
    % Outputs:
    %   model: trained model structure
    %   likelihood: variational lower bound
    %
    % Reference: Pattern Recognition and Machine Learning by Christopher M. Bishop (P.474)
    % Written by Mo Chen (sth4nth@gmail.com).
    %
    % Adapted for masked VB_GMM by Jordan Sorokin, 8/30/2017

    [d,n] = size(X);

    % initialize the prior hyperparameters
    if nargin < 4
        prior.alpha = 1; % dirichlet prior (mixing coefficients)
        prior.kappa = 1; % gaussian prior
        prior.mu = mean( X,2 ); % prior over means of gaussian clusters (i.e. one large cluster)
        prior.v = d+1; 
        prior.M = eye( d );   % M = inv(W) = covariance
    end
    prior.logW = -2 * sum( log( diag( chol( prior.M ) ) ) ); % inv(M) = conjugate prior for Sigma (multivariate covariance)

    % intialize our convergence parameters and presets
    tol = 1e-8;
    iter = 1;
    maxiter = 20;
    likelihood = -inf( 1,maxiter ); % the likelihood
    virtualX = compute_noise_ensemble( X',mask ); % computes the virutal distributions of the masked points in X
    model = initialize( X,m );
    
    % loop unitil convergence or maximum iteration
    while iter <= maxiter
        iter = iter+1;

        % run one step
        [model, lh] = VB_one_step( model );
        
        % check if best model
        if lh > max( likelihood(1:iter-1) ) && iter > 2
            bestModel = model;
            bestIter = iter;
        end
        
        % check for convergence
        likelihood(iter) = lh;
        dlh = abs( lh - likelihood(iter-1) ) / abs( lh );
        if dlh < tol
            break 
        end
    end

    % get final likelihood and labels
    %model = bestModel;
    %iter = bestIter;
    likelihood = likelihood(2:iter);
    [~,model.label(:)] = max( model.prob,[],2 );
    [~,~,model.label(:)] = unique( model.label );

    %% FUNCTIONS
    function [model,lh] = VB_one_step( model )
        % run one full step of the VB_EM algorithm

        % === E step ===
        model = expectation( virtualX.y',virtualX.eta',model );

        % === M step ===
        model = maximization( model );

        % === likelihoood ===
        lh = bound( model ) / n;
    end


    function model = initialize(X, m)
        n = size(X,2);
        if isstruct( m )  % initialize with a model
            model = m;
        elseif numel( m ) == 1  % random initialize k
            k = m;
            label = ceil( k*rand(1,n) ); % random labels up to k components
            model.prob = full( sparse( 1:n,label,1,n,k,n ) );
            model.label = label;
        elseif all( size( m )==[1,n] )  % initialize with labels
            label = m;
            k = max( label );
            model.prob = full( sparse( 1:n,label,1,n,k,n ) );
            model.label = label;
        else
            error('ERROR: init is not valid.');
        end
        model = maximization( model );
    end


    function model = expectation(X, eta, model)
        % computes the likelihood P( X | theta ) assuming the gaussian-wishart distribution

        %n = size( X,2 );
        [~,k] = size( model.mu );
        EQ = zeros( n,k );

        % loop over components 1:k
        for i = 1:k
            U = model.U(:,:,i); % U = chol( prior.M )
            invU = inv( U );
            Q = U' \ bsxfun( @minus,X,model.mu(:,i) ); 
            EQ(:,i) = d ./ (...
                model.kappa(i)... % dirichlet weight for this component
                + model.v(i) * dot( Q,Q,1 )... % X' * inv(Sigma) * X
                + sum( bsxfun( @times,eta,diag( invU * invU' ) ) )... % added normalization term
                ); % 10.64
        end

        ElogLambda = sum( psi( 0,0.5*bsxfun( @minus,model.v+1,(1:d)' ) ),1 ) + d*log( 2 ) + model.logW; % 10.65
        Elogpi = psi( 0,model.alpha ) - psi( 0,sum( model.alpha ) ); % 10.66 : 
        logRho = -0.5 * bsxfun( @minus,EQ,ElogLambda - d*log( 2*pi ) ); % 10.46
        logRho = bsxfun( @plus,logRho,Elogpi );   % 10.46
        model.logProb = bsxfun( @minus,logRho,logsumexp( logRho,2 ) ); % 10.49
        model.prob = exp( model.logProb );
        [~,model.label(:)] = max( model.prob,[],2 );
    end


    function model = maximization( model )
        % computes the prior distributions for each component
        % by updating model hyperparameters (i.e. parameters controlling prior 
        % distributions of the joint probability of latent parameters (mu,sigma) for 
        % each of the k components)

        nk = sum( model.prob,1 ); % 10.51
        k = numel( nk );
        alpha = prior.alpha + nk;
        kappa = prior.kappa + nk; % 10.60
        v = prior.v + nk; % 10.63
        Sigma = zeros( d,d );
        
        % masking logic
        p = sqrt( model.prob' );
        M = (mask' == 0); 
        L = ~M;
        
        % compute means (mean of unmasked points + virtual mean * unmasked points) 
        mu = bsxfun( @plus, prior.kappa*prior.mu, (virtualX.y'.*L + virtualX.mu'.*M)*model.prob); % do we need to divide by number in each cluster?
        mu = bsxfun( @rdivide,mu,kappa ); % 10.61: divides the mean by the sum of the probabilities...
        eta = (virtualX.eta'.*L + virtualX.sigma'.*M) * model.prob; % variances of masked features
        % mu = bsxfun( @plus,prior.kappa*prior.mu,virtualX.y'*model.prob );
        % mu = bsxfun( @times,mu,1./kappa ); % 10.61
        % eta = virtualX.eta' * model.prob;

        % update our current model
        model.U = zeros( d,d,k );
        model.logW = zeros( 1,k );
        model.alpha = alpha;
        model.kappa = kappa;
        model.mu = mu;
        model.v = v;
        mask_intersect = zeros( d,d );

        % loop over components
        for i = 1:k

            Xhat = bsxfun( @minus,virtualX.y',mu(:,i) );
            Xhat = bsxfun( @times,Xhat,p(i,:) );
            % Sigma = Xhat*Xhat';
            
            % loop over features for these points compute covariances
            for j = 1:d
               doubleM = M(j,:) & M; % masks only points that were masked for BOTH features i,j
               Sigma(j,:) =  Xhat(j,:) * (Xhat .* ~doubleM)'; 
               mask_intersect(:,j) = sum( doubleM,2 );
            end

            % now compute the modified covariance by taking our priors into account
            muhat = (prior.mu + virtualX.mu') - mu(:,i); % uses the prior mean AND the noise mean 
            Sigma_norm = (...
                 prior.M...
                 + Sigma... 
                 + prior.kappa * (mask_intersect .* (muhat*muhat'))...
                 + diag( eta(:,i) )...
                 );
            %  Sigma_norm = prior.M + Sigma + prior.kappa*(muhat*muhat');
            
            model.U(:,:,i) = chol( Sigma_norm );
            model.logW(i) = -2*sum( log( diag( model.U(:,:,i) ) ) );

            % Xhat = bsxfun( @minus,X,mu(:,i) ); % subtracts the mean of the ith component from unmasked points
            % Xhat = bsxfun( @times,Xhat,p(i,:) ); % multiplies by likelihood P( X | theta_k )
            % muhat = prior.mu - mu(:,i);
            % M = prior.M + Xhat*Xhat' + prior.kappa*(muhat*muhat');  % equivalent to 10.62 (prior cov. + new cov + hyperparameter kappa + cov of means  )
            % U(:,:,i) = chol( M );
            % logW(i) = -2*sum( log( diag( U(:,:,i) ) ) );      
        end
    end


    function likelihood = bound( model )
        k = size( model.prob,2 );
        Eqz = dot( model.prob(:),model.logProb(:) );
        logCalpha0 = gammaln( k*prior.alpha ) - k*gammaln( prior.alpha );
        Eppi = logCalpha0;
        logCalpha = gammaln( sum( model.alpha ) ) - sum( gammaln( model.alpha ) );
        Eqpi = logCalpha;
        Epmu = 0.5 * d * k * log( prior.kappa );
        Eqmu = 0.5 * d * sum( log( model.kappa ) );
        logB0 = -0.5 * prior.v * (prior.logW + d*log( 2 )) - logMvGamma( 0.5 * prior.v,d );
        EpLambda = k * logB0;
        logB =  -0.5 * model.v .* (model.logW + d*log( 2 )) - logMvGamma( 0.5 * model.v,d );
        EqLambda = sum( logB );
        EpX = -0.5 * d * n * log( 2*pi );
        likelihood = -Eqz+Eppi-Eqpi+Epmu-Eqmu+EpLambda-EqLambda+EpX;
    end
end 