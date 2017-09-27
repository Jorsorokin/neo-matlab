function mapplotMenu = st_plot_projections( handles )
    % function mapplotMenu = st_plot_projections( handles )
    %
    % plots the projected data onto the main axes of the sortTool GUI
    cla( handles.mapplot );
    delete( findobj( handles.figure1,'Type','uicontextmenu' ) ); 
    
    % plot a 2D scatter of the projections, and change the colors according
    % to the associated labels
    line( handles.mapplot,...
        handles.projection(:,handles.plotDims(1)),...
        handles.projection(:,handles.plotDims(2)),...
        'Marker','.','MarkerSize',30 );
    handles.mapplot.Children.CData = handles.plotcolor;
    set( gca,'color','k' );
    axis tight;

    % set up the right-click context menu capabilities 
    c = uicontextmenu;
    mapplotMenu = struct;
    handles.mapplot.UIContextMenu = c;
    mapplotMenu.rightClickDel = uimenu( 'Parent',c,'Label','Delete' );
    mapplotMenu.rightClickDel_label = uimenu( 'Parent',mapplotMenu.rightClickDel,'Label','label','Callback',@rightClick_delete );
    mapplotMenu.rightClickDel_point = uimenu( 'Parent',mapplotMenu.rightClickDel,'Label','point','Callback',@rightClick_delete );
    mapplotMenu.rightClickAssign = uimenu( 'Parent',c,'Label','Assign' );
    mapplotMenu.rightClickAssign_new = uimenu( 'Parent',mapplotMenu.rightClickAssign,'Label','new label','Callback',@rightclick_assign );
    mapplotMenu.rightClickAssign_existing = uimenu( 'Parent',mapplotMenu.rightClickAssign,'Label','existing label','Callback',@rightclick_assign );
end