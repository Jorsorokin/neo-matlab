function [labels,chisqStat,pval] = try_cluster_split( X,nDim )
    % [newLabels,chisqStat,pval] = try_cluster_split( X,nDim )
    %
    % Attempts to split the data in X into two clusters by modeling the
    % data as a chi-squared-distributed random variable with 
    % nDim-1 DOF. 
    %
    % X is first projected onto an nDim subspace using PCA,
    % which transforms the data to fit the assumptions of a 
    % chi-squared distribution (normally distributed,
    % independent variables). 
    %
    % Following the projection, the data is then modeled as a 
    % chi-squared variable with nDim-1 DOF as:
    %
    %       Z(i) = SUM_j{ (X(i,j) - X_hat(j))^2 }
    %               
    %   where X_hat(j) is the mean of dimension j across points.
    %
    % Under the null hypothesis, each Z(i) is a chi-square
    % random variable, and Z is distributed according to an 
    % nDim-1 chi-squared distribution. Then, a chi-square test statistic
    % can be implemented as:
    %
    %       chisq = SUM_i{ (O(i) - E(i))^2 / E(i) }
    %
    %   where O(i) and E(i) are observed and expected counts in the 
    %   empirical and theoretical CDFs of Z and an nDim-1 chi-squared
    %   distribution, respectively, for each bin i. 
    %
    % If the probability of observing the chisq statistic is less
    % than 0.025 (accounting for a two-tailed test), then the algorithm 
    % tries to split points in X via Spectral Clustering 
    % (more robust to odd shapes in the PCA space). The sub-clusters
    % are kept only if the sum of the mean squared error between the
    % points in each cluster to their own centers is less than the 
    % sum of the mean squared error of the points to the combined
    % cluster center. 
    % 
    % Inputs:
    %   X - an n x m matrix, with n = observations, m = dimensions
    %   
    %   nDim - the # of dimensions of the subspace to project onto
    %
    % Outputs:
    %   labels - a 1 x n vector of cluster labels
    %
    %   chisqStat - the chi squared statistic as described above
    %
    %   pval - the p-value of the chisqStat
    %
    % Written by Jordan Sorokin, 10/20/2017
    
    % check inputs, preallocate vars
    [n,m] = size( X );
    labels = ones( 1,n );
    errFunc = @(x,y)(sum( bsxfun( @minus,x,y ).^2,2 ));
    
    % (a) project the data using pca, to make the 
    % variables independent, normally distributed,
    % and maximally separated in the variance space
    [u,s] = svd( X','econ' );
    Y = X * (u(:,1:nDim) * s(1:nDim,1:nDim));

    % (b) sum the squared differences between the means of the
    % columns of Y and each row of Y (i.e take sum
    % of squared differences for each dimension)
    Y = Y ./ std( Y );          % each dimension is now unit variance...
    Yhat = mean( Y );           % mean of each dimension across points...
    Z = errFunc( Y,Yhat );      % Z is SUM_j{ (Y(i,j) - Yhat(j)) }
                                %   i.e. difference between each
                                %   point and the mean of points

    % (c) compare the distribution of Z with a chisquared
    % distribution of n DOF, where n = # of dimensions - 1.
    % If the points in Y come from one true multivariate gaussian 
    % distribution, then the resulting vector Z should be
    % approximately chisquare distributed with n DOF
    [chisqStat,pval] = check_chisquare_distribution( Z,nDim-1 );

    % (d) if the p-value of the chisquared statistic is small, perform
    % a split of the subclusters using spectral clustering, which
    % tends to be robust to oddly-shaped clusters that may occur
    % from splits of the PCA projection space
    if pval <= 0.005                

        % perform sepctral clusetering on the projected data
        [~,A] = adjacency_graph( Y,'kNN',10 );
        subClusts = SpectralClustering( A,2,2 );

        % compare new average waveforms to old average by
        % computing the mean squared error betwween points and centers
        Z_sub = zeros(1,2);
        for j = 1:2
            Y_sub = Y(subClusts == j,:); 
            %Y_sub = Y_sub ./ std( Y_sub );

            % compute the mean squared error of the
            % points-to-center of the new cluster
            Z_sub(1) = Z_sub(1) + mean( errFunc( Y_sub',mean( Y_sub )' ) );
            Z_sub(2) = Z_sub(2) + mean( errFunc( Y_sub',Yhat' ) ); 
        end
        
        % keep the clusters as separate clusters if the 
        % mean squared error of the sub clusters is smaller
        if Z_sub(1) < Z_sub(2)
            labels = subClusts;
        end        
    end
end