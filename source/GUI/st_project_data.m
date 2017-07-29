function [projection,W] = st_project_data( handles )
    % function [projection,W] = st_project_data( handles )
    %
    % Projects the data contained in "handles.data" onto a lower-dim subspace 
    % defined by "handles.projectMethod". "handles" refers to the structure
    % associted with the sortTool GUI
    %
    % returns the projected data and the projection matrix W, if applicable
    fprintf( 'Projecting data via %s\n',handles.model.projectMethod );
    pause( 0.5 );
    switch handles.model.projectMethod
        case 'ICA'
            W = fastica( handles.data',...
                'numOfIC',handles.nDim );
            projection = handles.W * handles.data;
        otherwise
            [projection,model] = compute_mapping( handles.data',...
                handles.model.projectMethod,handles.nDim );
            if isfield( model,'M' )
                W = model.M;
            else
                W = nan;
            end
    end

    % normalize the projections
    C = cov( projection );
    projection = (C^(-1/2) * projection')';   
end