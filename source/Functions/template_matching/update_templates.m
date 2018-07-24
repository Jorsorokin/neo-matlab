function templates = update_templates( templates,newTemplates )
    % templates = update_templates( templates,newTemplates )
    %
    % updates the template structure to reflect changes in templates or
    % completely new templates to be added
    
    remove = ismember( templates.labelMap,newTemplates.labelMap );
    templates.W(:,remove,:) = [];
    templates.M(:,remove) = []; 
    templates.prior(remove) = [];
    templates.w_norm(remove) = [];
    templates.A_norm.mu(remove,:) = [];
    templates.A_norm.sd(remove,:) = [];
    templates.labelMap(remove) = [];
    
    templates.W = [templates.W, newTemplates.W];
    templates.M = [templates.M, newTemplates.M];
    templates.prior = [templates.prior; newTemplates.prior*newTemplates.nPts / templates.nPts];
    templates.A_norm.mu = [templates.A_norm.mu; newTemplates.A_norm.mu];
    templates.A_norm.sd = [templates.A_norm.sd; newTemplates.A_norm.sd];
    templates.labelMap = [templates.labelMap; newTemplates.labelMap];
    templates.w_norm = [templates.w_norm; newTemplates.w_norm];
end