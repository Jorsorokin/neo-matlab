function fig = st_create_driftplot( handles )
    % fig = st_create_driftplot( handles )
    % creates an instance of the "driftplot" used by the "sortTool" GUI
    %
    % handles is a structure with variables referenced by the main GUI
    fig = figure( 'Visible','off','color','k','position',[700,80,800,250],'Visible','off' );
    set( fig,'MenuBar','none','Name','Spike drift' ); 
    
    ax = subplot( 'position',[.05 .1 .90 .5] ); hold on;
    set( ax,'tickdir','out','box','off','color','k',...
        'fontsize',8,'xcolor','w','ycolor','w',...
        'xlim',[min( handles.times ), max( handles.times )*max( handles.trials )] );        
    xlabel('time (s)');
    ylabel( 'amplitude' );
    
    ax = subplot( 'position',[.05 .65 .90 .3] ); hold on;
    set( ax,'color','k','ycolor','w','xcolor','k','xticklabel','','tickdir','out',...
        'xlim',[min( handles.times ), max( handles.times )*max( handles.trials )] );
    ylabel( 'num channels' );
end