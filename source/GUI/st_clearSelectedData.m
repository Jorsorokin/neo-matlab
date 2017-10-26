function handles = st_clearSelectedData( handles )
    % function handles = st_clearSelectedData( handles )
    %
    % clears the previously selected data (i.e. makes
    % "handles.selectedData = false), then updates plots
    if any( handles.selectedPoints )
        handles.mapplot.Children.SizeData(:) = 30;
        handles.selectedPoints(:) = false; % make all false

        % now update the plots
        st_clear_plots( handles );
    end
end