function ax = st_create_isiplot()
    % ax = st_create_isiplot()
    % creates an instance of the "isiplot" used by the "sortTool" GUI

    figure( 'Visible','off','color','k','position',[4,202,329,237] );
    set( gcf,'Name','Inter-spike intervals (ISI)' );
    ax = subplot( 1,1,1 ); hold on;
    set( ax,'xlim',[0 0.2] ); % 200 ms total
    set( ax,'tickdir','out','box','off','color','k','ylim',[0,1],...
        'ylimmode','auto','fontsize',8,'xcolor','w','ycolor','w' );
    set( ax,'LooseInSet',get( ax,'TightInSet' ) );
    xlabel( 'ISI (s)' );
end