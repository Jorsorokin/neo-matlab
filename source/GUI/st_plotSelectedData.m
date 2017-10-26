function st_plotSelectedData( handles )
    % function H = st_plotSelectedData( handles )
    %
    % plots the projected data points specified by "handles.selectedPoints"
    % onto the main axes of the sortTool GUI 
    if any( handles.selectedPoints )
        handles.mapplot.Children.SizeData(handles.selectedPoints) = 150;
    end
end