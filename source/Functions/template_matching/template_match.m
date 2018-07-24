function [ID,A] = template_match( X,W,minAmp,maxAmp )
    % [ID,A] = template_match( X,W,minAmp,maxAmp )
    %
    % performs greedy, iterative template matching based on the algorithm
    % proposed by Marre et al., 2012, J. Neurosci. 
    % The assumed model is that spikes are linear summation of parent
    % templates via the following equation:
    %
    %           Z_j = A_ij*W_j
    %   
    %           X_i = SUM{ Z_j + eta }
    %
    %   where A_ij are coefficients that scale the jth template
    %   waveform to match X. Most coefficients for
    %   the templates are actually set to 0, as any spike is unlikely to be
    %   a superposition of ALL possible templates in the case of many
    %   templates to begin with. Thus, the code limits the # of possible
    %   superpositions for any ith spike to 3 (see the nFitMax variable in
    %   the code, which can be changed).
    %
    % Inputs:
    %   X - an n x m matrix, with n = points, m = spikes
    %
    %   W - an n x k mean waveform template matrix, with k = # of templates
    %
    %   minAmp - minimum amplitude scaling factor for each template
    %
    %   maxAmp - maximum amplitude scaling factor for each template
    %
    % Outputs:
    %   ID - an m x 3 matrix of assigned labels. Note, one spike may be
    %        assigned more than one label, as in overlapping templats
    %   
    %   A - an m x 3 matrix of coefficients for W (see above)
    %
    % written by Jordan Sorokin, 6/7/2018
    
    % check sizes
    assert( size( X,1 ) == size( W,1 ),'Templates and data matrix X must have equal number of rows' );
    m = size( X,2 );
    k = size( W,2 ); 
    
    nFitMax = 3; % max 3 overlapping spikes for any given spike, OR 3 attempts at matching
    fullyMatchedBool = false( 1,m ); % keeps track of which spikes we'll iteratively match to templates
    templateCounter = zeros( 1,m,'uint8' );
    hasMatch = false( 1,m );
    tryTemplate = true( k,m ); % used for keeping track of (i,j) spike-template pairs
    A = zeros( m,nFitMax,class( X ) );
    ID = zeros( m,nFitMax,class( X ) );
    W_norm = zeros( 1,k,class( X ) );
    
    for j = 1:size( W,2 )
        W_norm(j) = norm( W(:,j,1) );
    end
    W_hat = W(:,:,1) ./ W_norm;
    
    Z = zeros( k,m );
    while any( templateCounter(~fullyMatchedBool) < nFitMax )
        
        % compute inner product with template waveform on un-matched spikes
        remainingPts = find( ~fullyMatchedBool );
        Z = W(:,:,1)' * X(:,remainingPts) .* tryTemplate(:,remainingPts);                                      
        
        % find the best matching template and update counters
        [a,bestTemplate] = max( Z );
        a = a ./ W_norm(bestTemplate); 
        inds = sub2ind( [k,m],bestTemplate,remainingPts );
        tryTemplate(inds) = false;
        templateCounter(remainingPts) = templateCounter(remainingPts) + 1;
        
        % mark invalid spikes as those outside of the amplitude bounds
        keep = (a > minAmp(bestTemplate,1)') & (a < maxAmp(bestTemplate,1)'); 
        keptPts = remainingPts(keep);
        X(:,keptPts) = X(:,keptPts) - a(keep).*W_hat(:,bestTemplate(keep)); % subtracts off the best-fitted mean waveform
        hasMatch(keptPts) = true;
        
        % loop over matched templates (faster than looping over matched
        % spikes) and update X
        uniqueTemplates = unique( bestTemplate(keep) );
        for i = uniqueTemplates
            idx = keep & (bestTemplate == i);
            matchedPts = remainingPts(idx);
            inds = sub2ind( [m,nFitMax],matchedPts,templateCounter(matchedPts) );
            A(inds) = a(idx);          
            ID(inds) = i;
            
            % subtract off the orthogonal components
            for j = 2:size( W,3 )
                temp = W(:,i,j)' * X(:,idx);
                X(:,idx) = X(:,idx) - temp.*(W(:,i,j)./norm( W(:,i,j) ));
            end
        end
            
        % update which spikes are fully matched 
        fullyMatchedBool(templateCounter == nFitMax) = true; % those that we've exhaustively searched
        fullyMatchedBool(find( hasMatch(remainingPts(~keep)) )) = true; % those that we've matched previously but couldn't match again
    end
    
    % finally, shift IDs over to left, filling up empty spaces
    for j = 1:2
        idx = ID(:,j) == 0 & ID(:,j+1) ~= 0;
        ID(idx,j) = ID(idx,j+1); ID(idx,j+1) = 0;
        A(idx,j) = A(idx,j+1); A(idx,j+1) = 0;
    end
    idx = ID(:,1) == 0 & ID(:,2) ~= 0;
    ID(idx,1) = ID(idx,2); ID(idx,2) = 0; 
    A(idx,1) = A(idx,2); A(idx,2) = 0;
end