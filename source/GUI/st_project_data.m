function [projection,mapping] = st_project_data( handles )
    % function [projection,mapping] = st_project_data( handles )
    %
    % Projects the data contained in "handles.data" onto a lower-dim subspace 
    % defined by "handles.projectMethod". "handles" refers to the structure
    % associted with the sortTool GUI
    %
    % returns the projected data and the projection matrix W, if applicable
    fprintf( 'Projecting data via %s\n',handles.R.projectMethod );
    
    switch handles.R.projectMethod
        case {'tSNE','LPP','Sammon','Isomap','NCA','LLE','HessianLLE','SNE','Autoencoder','MVU','fastMVU'}    
            [projection,~,mapping] = compute_spike_features( gather( permute( handles.data,[2,1,3] ) ),...
                handles.nDim,handles.R.projectMethod,handles.mask,handles.location,handles.concatChans );
        otherwise
            [projection,~,mapping] = compute_spike_features( permute( handles.data,[2,1,3] ),...
                handles.nDim,handles.R.projectMethod,handles.mask,handles.location,handles.concatChans );
    end
end