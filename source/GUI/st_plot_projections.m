function st_plot_projections( handles )
    % function st_plot_projections( handles )
    %
    % plots the projected data onto the main axes of the sortTool GUI
    cla( handles.mapplot );
    
    % plot a 2D scatter of the projections, and change the colors according
    % to the associated labels
    scatter( handles.mapplot,handles.projection(:,handles.plotDims(1)),...
        handles.projection(:,handles.plotDims(2)),30,'.');
    handles.mapplot.Children.CData = handles.plotcolor;
    set( gca,'color','k' );
    axis tight;  
end