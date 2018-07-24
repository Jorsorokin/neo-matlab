function fig = st_create_correlogramplot()
    % ax = st_create_correlogramplot()
    % creates an instance of the "correlogramplot" used by the "sortTool" GUI

    fig = figure( 'Visible','off','color','k','position',[600,300,500,400] );
    set( fig,'Name','Cross-correlograms' );
end