function st_plot_projections( handles )
    % function mapplotMenu = st_plot_projections( handles )
    %
    % plots the projected data onto the main axes of the sortTool GUI
    cla( handles.mapplot );
    
    % plot a 2D scatter of the projections, and change the colors according
    % to the associated labels
    scatter( handles.mapplot,...
        handles.projection(:,handles.plotDims(1)),...
        handles.projection(:,handles.plotDims(2)),...
        repmat(30,1,size( handles.projection,1 )),'.' );
    handles.mapplot.Children.CData = handles.plotcolor;
    handles.mapplot.Color = 'k';
    handles.mapplot.XColor = 'k';
    handles.mapplot.YColor = 'k';
    axis tight;
end