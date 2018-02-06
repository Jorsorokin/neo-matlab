function st_plot_projections( handles )
    % function mapplotMenu = st_plot_projections( handles )
    %
    % plots the projected data onto the main axes of the sortTool GUI
    delete( handles.mapplot.Children );
    
    scatterView = handles.projection * handles.plotDims;
    h = scatter( handles.mapplot,...
        scatterView(:,1),scatterView(:,2),...
        repmat(30,1,size( handles.projection,1 )),'.' );
    h.CData = handles.plotcolor;
    
    % set the bounds/color
    lims = [min( handles.projection(:) ),max( handles.projection(:) )];
    handles.mapplot.XLim = lims;
    handles.mapplot.YLim = lims;
    handles.mapplot.Color = 'k';
    handles.mapplot.XColor = 'k';
    handles.mapplot.YColor = 'k';
end