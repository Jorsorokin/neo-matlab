function model = get_sorting_model( X,labels,method,(mask) )
% model = get_sorting_model( X,labels,method )
%
% retrieve the model parameters used to sort the data in X 
% into distinct clusters. The clustering model depends on the actual
% sorting method used. 
%
% Inputs:
%	X - an n x d matrix of n, d-dimensional points
%	
%	labels - an n x 1 vector of integer labels, where labels(i) indicates 
%			 the cluster that the ith point is assigned to
%
%	method - a string specifying the cluster method used. Valid methods are:
%			'EM-GMM'	:	expectation-maximization for gaussian mixture model
%			'Km'		:	k-means
%			'VB'		:	variational bayes for gaussian mixture model
%			'DBSCAN'	:	density-based spatial clustering 
%			'Spectral'	:	spectral clustering using the eigenmap of the laplacian matrix
%
% Outputs:
%	model - a structure containing the necessary parameters for future clustering using 
%			the same sorting method 
%
% By JMS, 9/15/2017

model = struct;

% check the method
switch method
	case {'EM-GMM','VB','Km'}
		[model.mu,model.Sigma] = get_cluster_description( X,labels );
		[model.probabilities,model.w] = get_cluster_probabilities( X,labels,model.mu,model.Sigma );

	case 'DBSCAN'
		fprintf( 'DBSCAN model creation under development\n' );

	case 'Spectral'
		fprintf( 'Spectral clustering model creation under development\n' );
end


end
