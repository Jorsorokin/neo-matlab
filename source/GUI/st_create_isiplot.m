function fig = st_create_isiplot()
    % fig = st_create_isiplot()
    % creates an instance of the "isiplot" used by the "sortTool" GUI

    fig = figure( 'Visible','off','color','k','position',[4,202,329,237],'name','Inter-spike interval (ISI)' );
    ax = subplot( 1,1,1 );
    set( ax,'tickdir','out','box','off','color','k','xlim',[0.1 1000],'ylim',[0 1],...
        'ylimmode','auto','fontsize',8,'xcolor','w','ycolor','w','xscale','log','NextPlot','add' );
    set( ax,'LooseInSet',get( ax,'TightInSet' ) );
    xlabel( 'ISI (ms)' );
    ylabel( 'P(x)' );
end