function [labels,labelMap] = refine_clusters( X,labels,varargin )
    % [labels,labelMap] = refine_clusters( X,labels,(mask,merge,nDim) )
    %
    % Refines the clusters specified by "labels" by first attempting to 
    % split each subset of data in X belonging to a unique cluster,
    % then attempting to merge clusters together. This is a wrapper 
    % that calls "try_cluster_split.m". 
    % Please see "try_cluster_split.m" for details on the algorithm and implementation
    %
    % Inputs:
    %   X - n x m x c matrix, with n = pts, m = observations, c = channels
    %
    %   labels - 1 x m vector of cluster labels
    %
    %   (mask) - c x m masking matrix, specifying which channels
    %            were "informative" for each observation 1:m
    %
    %   (merge) - logical flag. If true, will attempt to merge pairs of
    %             clusters together. Default = false
    %
    %   (nDim) - the number of dimensions to project each cluster onto via
    %            PCA. Default = 3;
    %
    % Outputs:
    %   labels - 1 x n updated label vector
    %
    %   labelMap - 1 x K vector, mapping each new label to the previous 
    %              label. For instance, labelMap(:,4) = [20,21] means
    %              the cluster assigned label = 4 was split into two new 
    %              clusters with labels 20 and 21
    %
    % Written by Jordan Sorokin, 10/20/2017
    
    
    % check inputs
    if nargin > 2 && ~isempty( varargin{1} ) && ~any( isnan( varargin{1} ) )
        mask = varargin{1};
    else
        mask = [];
    end
    if nargin > 3 && ~isempty( varargin{2} )
        merge = varargin{2};
    else
        merge = false;
    end
    if nargin > 4 && ~isempty( varargin{3} )
        nDim = varargin{3};
    else
        nDim = 3;
    end
    
    allLabels = unique( labels(labels>0) );
    maxLabel = max( allLabels );
    labelMap = zeros( 2,maxLabel );
    split = false( 1,maxLabel );
    
    %% Part 1: loop over clusters and try to split each
    fprintf( 'Attempting to split clusters...\n' );
    for k = allLabels
        try
            pts = find( labels == k );
            
            % if mask is not nan, only use channels that are not masked
            if ~isempty( mask )
                keepChans = any( mask(:,pts)' > 0 ); % only the most informative channels
                X_k = concatenateSpikes( X(:,pts,keepChans) );
            else
                X_k = concatenateSpikes( X(:,pts,:) );
            end

            % try splitting the cluster. Update the label vector
            % if splitting was successful
            newID = try_cluster_split( X_k',nDim );
            if any( newID == 2 )
                split(k) = true;
                for j = 1:2
                    maxLabel = maxLabel + 1;
                    labelMap(j,k) = maxLabel;
                    labels(pts(newID==j)) = maxLabel;
                end
            end            
        catch
            fprintf( 'error trying to split cluster %i\n',allLabels(k) );
            disp( lasterr );
        end
    end
    
    %% Part 2: try to merge clusters via pairwise comparisons
    if merge
        fprintf( 'Attempting to merge clusters...\n' );
        triedMerge = false( maxLabel );
        merged = false( 1,maxLabel );
        allLabels = unique( labels(labels>0) );
        for j = find( split )
            triedMerge(labelMap(1,j),labelMap(2,j)) = true;
            triedMerge(labelMap(2,j),labelMap(1,j)) = true;
        end

        for k = allLabels 
            % skip this cluster if it was merged
            if merged(k)
                continue
            end

            pts_k = find( labels == k );

            for j = allLabels(allLabels ~= k)
                if triedMerge(j,k) || merged(j)
                    continue
                end

                % pull out the points for these two clusters
                triedMerge(k,j) = true;
                pts_j = find( labels == j );

                % if mask exists, only use unmasked channels
                if ~isempty( mask )
                    X_kj = concatenateSpikes( X(:,[pts_k,pts_j],any( mask(:,[pts_k,pts_j])' > 0 )) );
                else
                    X_kj = concatenateSpikes( X(:,[pts_k,pts_j],:) );
                end

                % Try merging the clusters and update the labelMap if we do.
                % Although we're calling "try_cluster_split", we only extract
                % the chisqStat and pval from the merged chi-square test, and 
                % can use this to see if the merged data follow a chi-squared
                % distribution around their mean
                [~,~,p] = try_cluster_split( X_kj,nDim );
                if p >= 0.1
                    merged([k,j]) = true;

                    % find the previous if these clusters were split previously
                    switch min( k,j )
                        case j
                            labels(pts_k) = j;
                            prevK = (labelMap == k);                
                            if ~any( prevK )
                                labelMap(1,k) = j;
                            else
                                labelMap(prevK) = j;
                            end
                        otherwise
                            labels(pts_j) = k;
                            prevJ = (labelMap == j);                
                            if ~any( prevJ )
                                labelMap(1,j) = k;
                            else
                                labelMap(prevJ) = k;
                            end
                    end
                end
            end
        end
    end
end
        
        