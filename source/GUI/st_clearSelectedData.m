function handles = st_clearSelectedData( handles )
    % function handles = st_clearSelectedData( handles )
    %
    % clears the previously selected data (i.e. makes
    % "handles.selectedData = false), then updates plots
    axes( handles.mapplot );
    handles.selectedPoints(:) = false; % make all false
    if isfield( handles,'selectedPointPlot' ) 
        if ~isa( handles.selectedPointPlot,'double' ) && isvalid( handles.selectedPointPlot )
        
            % delete the previously selected points
            delete( handles.selectedPointPlot );

            % now update the plots
            st_clear_plots( handles );
        end
    end
end