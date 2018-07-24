function fig = st_create_rasterplot( handles )
    % fig = st_create_rasterplot( handles )
    % creates an instance of the "rasterplot" used by the "sortTool" GUI
    %
    % handles is a structure with variables referenced by the main GUI
    fig = figure( 'Visible','off','color','k','position',[260,80,1000,200] );
    set( fig,'MenuBar','none','Name','Spike rasters' ); 
    
    ax = subplot( 'position',[.05 .1 .95 .65] ); hold on;
    set( ax,'tickdir','out','box','off','color','k',...
        'fontsize',8,'xcolor','w','ycolor','w','ytick',1:max( handles.trials ),...
        'ylim',[min(handles.trials)-1, max( handles.trials)],...
        'xlim',[min(handles.times),max(handles.times)] );        
    xlabel('time (s)');
    ylabel('trials');
    
    ax = subplot( 'position',[.05 .8 .95 .15] ); hold on;
    set( ax,'color','k','ycolor','k','xcolor','k',...
        'yticklabel','','xticklabel','',...
        'xlim',[min(handles.times),max(handles.times)] );
end