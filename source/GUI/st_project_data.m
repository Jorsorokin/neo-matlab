function [projection,mapping] = st_project_data( handles,varargin )
    % function [projection,mapping] = st_project_data( handles )
    %
    % Projects the data contained in "handles.data" onto a lower-dim subspace 
    % defined by "handles.projectMethod". "handles" refers to the structure
    % associted with the sortTool GUI
    %
    % returns the projected data and the projection matrix W, if applicable
    fprintf( 'Projecting data via %s\n',handles.model.projectMethod );
    pause( 0.5 );
    
    % check if we should concatenate or not (i.e. for masked-EM, no
    % concatenation)
    if nargin > 1 && ~isempty( varargin{1} )
        concatenate = varargin{1};
    else
        concatenate = true;
    end
    
    switch handles.model.projectMethod
        case {'tSNE','LPP','Sammon','Isomap','NCA','LLE','HessianLLE','SNE','Autoencoder','MVU','fastMVU'}    
            [projection,~,mapping] = compute_spike_features( gather( permute( handles.data,[2,1,3] ) ),...
                handles.nDim,handles.model.projectMethod,handles.mask,handles.location,concatenate );
        otherwise
            [projection,~,mapping] = compute_spike_features( permute( handles.data,[2,1,3] ),...
                handles.nDim,handles.model.projectMethod,handles.mask,handles.location,concatenate );
    end
end