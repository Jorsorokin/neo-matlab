function [chisqStat,pval,V,Y] = chi2normal( X )
    % [chisqStat,pval,V,Y] = chi2normal( X )
    %
    % determines whether X is multivariate-normally distributed using the
    % chisquared goodness of fit test. This is a wrapper that first
    % whitens the columns of X (so that they are independent and have equal
    % variances) and then performs a chi-square test on the result. Under
    % the null, the columns in X are all gaussian distributed and indeed
    % the distribution of X is unimodal, multi-variate normal. 
    %
    % Thus, this is NOT a test that the columns of X are dependent, since
    % they are first whitened. Insead, it tests whether the multi-variate
    % data shows non-gaussianity, such as a bi-modal distribution. 
    %
    % Inputs:
    %   X - n x m matrix, with n = observations, m = variables
    %
    % Outputs:
    %   chisqStat - the chisquared test statistic
    %
    %   pval - the p-value of having observed a particular goodness of fit 
    %          given the null (X is multivariate normal)
    %
    %   V - cramer's V effect size (sqrt( chisqStat / n )). Helps decide
    %       how large of an "effect" (in this case, how different from a
    %       chi-squared distribution X is). Values > 0.3 are considered
    %       moderate to large
    %
    %   Y - the whitened data
    %           Y = X * (cov( X )^-0.5)' (zca whitening)
    %
    % Written by Jordan Sorokin, 7/18/18
    
    [n,m] = size( X );
    
    % (a) standardize the dimensions of X, a necessary whitening 
    % step for the chisquare distribution which assumes independent,
    % normally distributed data
    %W = (cov( X )^-0.5 + eye( m )*1e-8)';
    %Y = X*W;
    Y = zscore( X );
    
    % (b) sum the squared differences between the means of the
    % columns of Y and each row of Y (i.e take sum
    % of squared differences for each dimension)
    Z = sum( (Y - mean( Y )).^2,2 );

    % (c) compare the distribution of Z with a chisquared
    % distribution of n DOF, where n = # of columns. 
    % If the points in Y come from one true multivariate gaussian 
    % distribution, then the resulting vector Z should be
    % approximately chisquare distributed with n DOF
    [chisqStat,pval] = check_chisquare_distribution( Z,m );
    
    % (d) compute cramers V 
    V = sqrt( chisqStat / n );
end