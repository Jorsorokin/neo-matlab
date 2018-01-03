function ax = st_create_rasterplot( handles )
    % ax = st_create_rasterplot( handles )
    % creates an instance of the "rasterplot" used by the "sortTool" GUI
    %
    % handles is a structure with variables referenced by the main GUI
    figure( 'Visible','off','color','k','position',[260,80,1000,280] );
    set( gcf,'MenuBar','none','Name','Spike rasters' ); 
    ax = subplot( 1,1,1 ); hold on;
    set( ax,'tickdir','out','box','off','color','k',...
        'fontsize',8,'xcolor','w','ycolor','w',...
        'ylim',[min(handles.trials)-1, max( handles.trials)+1],...
        'xlim',[min(handles.times),max(handles.times)] );
    set( ax,'LooseInSet',get( ax,'TightInSet' ),'position',[.05 .15 .95 .8] );
    xlabel( 'time (s)' );
    ylabel( 'trials' )
end