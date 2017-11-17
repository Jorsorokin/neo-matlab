function quality = measure_cluster_quality( X,labels,varargin )
	% quality = measure_cluster_quality( X,labels,(method,distMeasure) )
	%
	% Evaluates the clustering of the points in X. Interpretation of 
	% cluster quality depends on the method used (see specific methods
	% for interpretation of the quality)
	% 
	%
	% Inputs:
	%	X - npts x ndim matrix of data points 
	%
	%	labels - npts-length vector of cluster IDs of the points in X
	%
	%	(method) - string specifying the evaluation metric
	%				'Silhouette' (default)
	%				'Calinski-Harabasz'
	%				'Davies-Bouldin'
	%	
	%	(distMeasure) - string or user-defined function specifying the distance measure to use.
	%					See "pdist2" for valid distance measurements. Default = 'euclidean'
	%				
	% Outputs:
	%	quality - 1 x K vector of cluster quality, for each cluster K.
	%			
	% Written by Jordan Sorokin
	% 11/16/2017

	% check inputs
	[n,m] = size( X );
	if n ~= length( labels )
		error( '# of points in X does not match # of points in the label vector' );
	end

	if nargin > 2 && ~isempty( varargin{1} )
		method = varargin{1};
	else
		method = 'Silhouette';
	end

	if nargin > 3 && ~isempty( varargin{2} )
		distMeasure = varargin{2};
	else
		distMeasure = 'euclidean';
	end

	clusts = unique( labels );
	K = numel( clusts );

	% switch over the evaluation metric
	switch method
		case 'Silhouette'

			quality = zeros( 1,K );
			
			% compute the avg distance between point i and other pts
			for k = 1:K
				pts = (labels == clusts(k));
				a = mean( pdist2( X(pts,:),X(pts,:),distMeasure ),2 );
				b = inf( nnz( pts ),1 );

				for j = clusts(~ismember( clusts,k ))
					pts2 = (labels == j);
					d = mean( pdist2( X(pts,:),X(pts2,:),distMeasure ),2 );
					b = min( b,d );
				end

				% store silhouette index for this cluster
				quality(k) = mean( (b - a) ./ max( a,b ) );
			end

		case 'Calinski-Harabasz'

			c = mean( X ); % overall centroid
			W = zeros( 1,K ); % within-cluster distance
			B = zeros( 1,K ); % between-group dispersion

			% compute within-cluster dispersions
			for k = 1:K
				pts = (labels == clusts(k));
				c_k = mean( X(pts,:) );
				W(k) = sum( pdist2( X(pts,:),c_k,distMeasure ) );
				B(k) = pdist2( c_k,c,distMeasure );
			end

			% store dispersion ratios
			quality = (B - W) ./ W .* (n - K / (K - 1));

		case 'Davies-Bouldin'
			
			S = zeros( 1,K );
            A = zeros( K,m,class( X ) );
            
			for k = 1:K
				pts = (labels == clusts(k));
				A(k,:) = mean( X(pts,:) );
				S(k) = 1/nnz( pts ) * sum( pdist2( X(pts,:),A(k,:),distMeasure ) );
			end

			M = pdist2( A,A,distMeasure );
			R = (S + S') ./ M;
			R = R .* (ones( K,K )+ -1e6*eye( K ));
			
			quality = max( R );
	end
end







