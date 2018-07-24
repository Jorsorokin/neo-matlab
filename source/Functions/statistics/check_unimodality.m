function [D,mu,sigma,cutloc,poly] = check_unimodality( x,varargin )
    % D = check_unimodality( x ) determines whether the
    % input vector x has a uni-modal distribution by taking sequential cuts
    % along the empirical cdf of x and comparing the means of the lower / upper
    % portions as:
    %   
    %       D(i) = ( mu(i,1) - mu(i,2) ) / ( sigma(i,1) + sigma(i,2) )
    %
    %   where mu and sigma are the means and SDs of the two populations
    %   following the ith cut
    %
    % If x is unimodal, D will be convex with a local minimum at its mean. If x
    % is not unimodal, D will be concave, with a peak at the point that
    % maximally separates the two true distributions that comprise x. If there
    % is more than one peak, this indicates x may be n-modal, where n is
    % determined by the number of peaks in D.
    %
    % [D,mu,sigma,cutloc] = check_unimodality( x ) also returns the means and SDs 
    % of the populations at the various cuts, and the bin location of each cut
    %
    % [D,mu,sigma,cutloc] = check_unimodality( x,bins ) allows the user to specify the
    % bins for the histogram prior to taking the cumulative distribution
    % (default is to use 51 evenly spaced bins between the range of x)
    %
    % [D,mu,sigma,cutloc] = check_unimodality( ...,cuts ) also allows the user to
    % specify the location of the cuts along the x-axis of the distribution
    % (default is to cut along 5% increments in the total data distribution).
    % "cuts" should be defined as percentage of the total cdf
    %
    % [...,poly] = check_unimodality( ... ) also returns a 2nd degree
    % polynomial fit structure to "D"
    %
    % by Jordan Sorokin, 4/4/18

    % check inputs
    if nargin > 1 && ~isempty( varargin{1} )
        bins = varargin{1};
    else
        bins = [];
    end
    
    if nargin > 2 && ~isempty( varargin{2} )
        cuts = varargin{2};
        nCuts = numel( cuts );
    else
        nCuts = 19;
        cuts = linspace( 0.05,0.95,nCuts );
    end
    
    % create the pdf & cdf
    if isempty( bins )
        [P,bins] = histcounts( x );
    else
        P = histcounts( x,bins );
    end
    
    C = cumsum( P ) / sum( P );

    % cut C along various points and compute D(i) for the ith cut
    mu = zeros( nCuts,2 );
    sigma = zeros( nCuts,2 );
    cutloc = zeros( nCuts,1 );

    for i = 1:nCuts
        idx = find( C > cuts(i),1 );
        cutloc(i) = bins(idx);
        
        xsub1 = x(x < cutloc(i)); % lower population
        xsub2 = x(x > cutloc(i)); % upper population

        mu(i,1) = mean( xsub1 );
        mu(i,2) = mean( xsub2 );
        sigma(i,1) = std( xsub1 );
        sigma(i,2) = std( xsub2 );
    end

    D = ( mu(:,2) - mu(:,1) ) ./ sum( sigma,2 );
    
    % now fit a polynomial to D to smooth out noise 
    loc = linspace( cutloc(1),cutloc(end),100 );
    [coef,S,m] = polyfit( cutloc,D,2 );
    [fit,err] = polyval( coef,loc,S,m );

    % determine if convex or concave, then find trough / peak
    dfit = diff( fit );
    if dfit(1) < 0 && dfit(end) > 0
        [pk,pt] = min( fit );
        unimodal = true;
    else
        [pk,pt] = max( fit );
        unimodal = false;
    end
    poly = struct( 'fitX',loc,'fitY',fit,'fitErr',err,'fitCoef',coef,'fitStruct',S,'peakVal',pk,'peakLoc',loc(pt),'unimodal',unimodal );
end