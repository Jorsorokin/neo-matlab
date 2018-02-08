function [labels,model] = predict_clusters( X,method,model )
% function [labels,model] = predict_clusters( X,method,model )
%
% assign cluster labels to each data point in X given a previous
% sorting model / method
%
% Inputs:
%   X       - n x m data matrix, with n = observations, m = variables
%
%   method  - string specifying the clustering method. Valid strings are:
%               'EM-GMM' - expectation maximization gaussian mixture model
%               'VB' - variational bayes gaussian mixture model
%               'Km' - kmeans
%               'DBSCAN' - density based spatial clustering
%               'Spectral' - spectral clustering
%
%   model   - sorting model structure with the necessary fields for the given method
%
% Outputs:
%   labels  - 1 x n vector of cluster assignments
%
%   model   - (potentially) updated structure given the new data points in X
%
% Written by Jordan Sorokin, 9/17/2017

% check inputs
if ~isstruct( model )
    error( 'model is not a structure. Please see "sort_clusters.m" for details' );
end

% switch of the type of sorting routine, and perform cluster assignment    
fprintf( 'Predicting cluster labels via %s\n',method );
pause(0.2);

switch method
    case {'EM-GMM','VB'}
        if strcmp( method,'EM-GMM' )
            labels = mixGaussPred( X',model );
        otherwise
            labels = mixGaussVbPred( X',model );
        end

    case 'Km'
        % get the closest centroid to each data point
        labels = nn_batch( X,model.mu,1 );

    case 'Spectral'
        W = adjacency_graph( X,model.neighborType,model.nNeighbors );
        D = diag( sum( W,2 ) ); 
        labels = wMSC( W,D,model.sigma,model.K )

    case 'DBSCAN'
        fprintf( 'DBSCAN cluster prediction under development...\n' );
        labels = nan;
        
    case 'HDBSCAN'
        labels = model.predict( X );

    otherwise
        fprintf( 'Unable to predict cluster assignment via %s clustering.\n',method );
        labels = nan;
end

end
