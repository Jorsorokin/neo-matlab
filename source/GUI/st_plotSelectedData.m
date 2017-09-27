function H = st_plotSelectedData( handles )
    % function H = st_plotSelectedData( handles )
    %
    % plots the projected data points specified by "handles.selectedPoints"
    % onto the main axes of the sortTool GUI as a "*"
    if any( handles.selectedPoints )
        
        % delete previous points
        if isfield( handles,'selectedPointPlot' ) 
            if ~isa( handles.selectedPointPlot,'double' ) && isvalid( handles.selectedPointPlot )
                delete( handles.selectedPointPlot );
            end
        end
        
        % now create the new points
        axes( handles.mapplot ); hold on;
        colors = handles.mapplot.Children.CData(handles.selectedPoints,:);
        plotPts = handles.projection(handles.selectedPoints,handles.plotDims);
        H = scatter( plotPts(:,1),plotPts(:,2),50,'*' );
        H.CData = colors;
    else
        H = nan;
    end
end