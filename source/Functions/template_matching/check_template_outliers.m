function outliers = check_template_outliers( templates )
    % outliers = check_template_outliers( templates )
    %
    % checks if any template is an "outlier" by comparing the variance of
    % the amplitude fits to the distribution of variances. Templates with
    % high variance fits are assumed to have been created by possibly
    % multiple putative neurons, each of which fits poorly with the
    % template
    
    v = templates.A_norm.sd(:,1);
    outliers = find( v > (median( v ) + 2*mad( v )) )';
end