function st_clearSelectedData( handles )
    % function st_clearSelectedData( handles )
    %
    % clears the previously selected data (i.e. makes
    % "handles.selectedData = false), then updates plots
    if any( handles.selectedPoints )
        handles.mapplot.Children.SizeData(handles.selectedPoints) = 30;
        st_clear_plots( handles );
    end
end