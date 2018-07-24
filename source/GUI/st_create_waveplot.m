function fig = st_create_waveplot( handles )
    % fig = st_create_waveplot( handles )
    % creates an instance of the "waveformplot" used for the "sortTool" GUI
    %
    % handles is a structure with variables referenced by the main GUI
    fig = figure( 'Visible','off','position',[4,500,500,400],'color','k','Name','Raw waveforms' );
    if handles.availableData
        yrange = gather( [min( min( handles.data(:) ) ),max( max( handles.data(:) ) )] );
        xrange = gather( [1, size( handles.data,1 )] );
    else
        yrange = [-500 200];
        xrange = [0 1];
    end
    nchan = size( handles.data,3 );
    cols = ceil( sqrt( nchan ) );
    rows = ceil( nchan / cols );
    marg_h = [0.03 0.03]; marg_w = [0.04 0.04]; gap = [0.01 0.01];
    axh = (1-sum(marg_h)-(rows-1)*gap(1))/rows; % x positions
    axw = (1-sum(marg_w)-(cols-1)*gap(2))/cols; % y positions
    py = 1-marg_h(2)-axh; 
    
    for ih = 1:rows
        px = marg_w(1);
        for ix = 1:cols
            ax = axes( 'Units','normalized', ...
                'Position',[px py axw axh], ...
                'XTickLabel','',...
                'YTickLabel','',...
                'NextPlot','add' );
            
            set( ax,'ylim',yrange,'xlim',xrange,...
                'fontsize',8,'ycolor','k',...
                'xcolor','k','color','k' );
            
            px = px+axw+gap(2);
        end
        py = py-axh-gap(1);
    end
    linkaxes( get( gcf,'Children'),'xy' );
end