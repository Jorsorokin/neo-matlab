function fig = st_create_dvdtplot( handles )
    % fig = st_create_driftplot( handles )
    % creates an instance of the "dvdtplot" used by the "sortTool" GUI
    %
    % handles is a structure with variables referenced by the main GUI
    fig = figure( 'Visible','off','color','k','position',[800,500,250,250],'Visible','off' );
    set( fig,'MenuBar','none','Name','Phase plot' ); 
    
    ax = subplot( 'position',[.1 .1 .85 .85] );
    set( ax,'tickdir','out','box','off','color','k',...
        'ycolor','w','xcolor','w','NextPlot','add' );        
    xlabel('V');
    ylabel( 'dVdT' );
end