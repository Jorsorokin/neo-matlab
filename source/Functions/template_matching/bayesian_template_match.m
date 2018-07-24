function [labels,P] = bayesian_template_match( X,Y,model,mask,alpha,beta )
    % [labels,P] = bayesian_template_match( X,Y,model,mask )
    %
    % performs bayesian template matching on the data matrix X given the template model
    % created by "create_spike_templates.m"
    %
    % Templates Y_i are fit to each data point X_j in a probabilistic framework:
    %
    %                 P(Y_i | X_j) = P(X_j | Y_i) * P(Y_i)
    %       
    %       - P(X_j | Y_i) = P(theta | Y_i) * (e^-(|| (X_j(t) - Y_i(t)) / 2*sigma(t)^2 ||^2)
    %       - P(theta | Y_i) = P(peak time,channel location,half width | Y_i)
    %       - P(Y_i) = prior probability of cluster i
    %
    % Rather than formulating the above as a multi-variate gaussian
    % distribution, the sample of points used to create the templates in
    % the first place is used to estimate the histogram of possible values
    % for the various parameters in "theta". 
    %
    % Because some new points may be  beyond the limits of the histogram, 
    % we also use the second term in P(X_j | Y_i) - the SSE between X_j and
    % Y_i - as another factor to avoid 0 proabilities. "sigma" is defined
    % during template creation as the variance for each time point (t)
    % across all points used to create the ith template
    %
    % By Jordan Sorokin, 6/4/2018
    
    % check inputs
    [n,m,c] = size( X );
    [n2,c2,K] = size( Y );
    assert( all( [n,c] == [n2,c2] ),...
        'Templates and data matrix X must have equal # of points and channels' );
    
    Y = reshape( Y,n*c,K );
    sigma = reshape( model.sigma,n*c,K );
    P_y = log( model.prior );
    
    L = full( (1:c) * mask ./ sum( mask,1 ) );
    C = full( sum( mask > 0,1 ) );
    [~,bestChan] = max( mask,[],1 );
    H = zeros( 1,m );
    A = zeros( 1,m );
    T = zeros( 1,m );
    for ch = 1:c
        idx = (bestChan == ch);   
        if nnz( idx ) == 0
            continue
        end
        [A(idx),T(idx)] = min( X(:,idx,ch) );
        H(idx) = halfWidth( X(:,idx,ch),A(idx),1 );
    end
    
    % pre-allocate vars
    P = ones( m,K,class( X ) );
    parameters = fields( model )';
    X = concatenateSpikes( X );
    MU = [];
    VAL = [];
    covar = model.covar;
    
    % loop over templates
    for param = parameters
        switch param{1}
            case 'amplitude'
                val = A';
            case 'peaktime'
                val = T';
            case 'halfwidth'
                val = H';
            case 'location'
                val = L';
            %case 'nchans'
            %    val = C';
            otherwise
                continue
        end 
        
        evalc( sprintf('mu_i = model.%s.mu',param{1}) );
        evalc( sprintf('sigma = model.%s.sd',param{1}) );
        evalc( sprintf('xvec = model.%s.x',param{1}) );
        evalc( sprintf('probability = model.%s.p',param{1}) );
        
        [~,closestPt] = min( compute_pairwise_dist( xvec',val ) ); 
        
        for i = 1:K
            %P(:,i) = P(:,i) - (val - mu(i)).^2 / (2*sigma(i)^2 + eps);
            P(:,i) = P(:,i) + log( probability(i,closestPt) + eps )';
        end
    end
    
    % now add L2-norm between X and all Y
    P_xy = -compute_pairwise_dist( X',Y' );
    
    % find cluster that maximizes the probability
    P = beta*(alpha*P + (1-alpha)*P_xy) + (1-beta)*P_y';
    [~,labels] = max( P,[],2 );
end