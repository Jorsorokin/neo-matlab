function [model,K,labels,posterior,bic] = fit_optimum_gmm( X,maxK,varargin )
    % [model,K,labels,posterior,bic] = fit_optimum_gmm( X,maxK,(lambda,replicates,covtype,sharedcov,outlierpdf) )
    %
    % sort using GMM with up to maxK-components, then find best model via BIC
    % minimization. Specify optional arguments using the name-value pair
    %
    % Inputs:
    %   X - an N x D matrix, with N = observationd, D = dimensions/variables
    %
    %   maxK - integer specifying max # of clusters to search for
    %
    %   (lambda) - decimal specifying regularization parameter to add to
    %              diagonal of the covariance matrix (default = 0.01)
    %
    %   (replicates) - # of replicate fits to perform to avoid local optima
    %                  (default = 1)
    %
    %   (covtype) - the covariance type for the data (see fitgmdist.m)
    %               (default = 'full')
    %
    %   (sharedcov) - boolean specifying shared or not shared covariance 
    %                 (default = false)
    %
    %   (outlierpdf) - the minimum pdf of the points evaluated via the best 
    %                  found model below which points are considered outliers
    %                  (default = 1e-7)
    %
    % Outputs:
    %   model - the best selected model
    %
    %   K - the best-estimated # of components in the data
    %
    %   labels - the result of hard clustering followed by outlier removal
    %
    %   posterior - an N x K matrix of posterior probabilities of the points
    %
    %   bic - a 1 x maxK vector of BIC values (best model is the minimum)
    %
    % Written by Jordan Sorokin
    % 7/16/18

    % get inputs
    p = check_inputs( varargin );

    % loop over gmm models and fit data with each
    gmmModels = cell( 1,maxK );
    bic = inf( 1,maxK );
    for k = 1:maxK
        gmmModels{k} = fitgmdist( X,k,...
                                  'Regularization',p.regularization,...
                                  'Replicates',p.replicates,...
                                  'CovarianceType',p.covtype,...
                                  'SharedCovariance',p.sharedcov );
        bic(k) = gmmModels{k}.BIC;
    end

    % best model minimizes the BIC
    [~,K] = min( bic );
    model = gmmModels{K};

    % detect outliers using the posteriors to clean up the
    % clusters / remove noise clusters
    [labels,~,posterior,logpdf] = model.cluster( X );
    outliers = exp( logpdf ) <= p.outlierpdf;
    labels(outliers) = 0;
    
    %% HELPER
    function p = check_inputs( inputs )
        names = {'regularization','replicates','covtype','sharedcov','outlierpdf'};
        defaults = {0.01,1,'full',false,1e-7};
        
        p = inputParser();
        for i = 1:numel( names )
            p.addParameter( names{i},defaults{i} );
        end
        
        p.parse( inputs{:} );
        p = p.Results;
    end       
end