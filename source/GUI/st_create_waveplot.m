function ax = st_create_waveplot( handles )
    % ax = st_create_waveplot( handles )
    % creates an instance of the "waveformplot" used for the "sortTool" GUI
    %
    % handles is a structure with variables referenced by the main GUI
    figure( 'Visible','off','position',[4,530,340,250],'color','k' );
    set(gcf,'Name','Raw waveforms');
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
    
    rowCounter = 0;
    ii = 0;
    for ih = 1:rows
        px = marg_w(1);

        for ix = 1:cols
            ii = ii+1;
            ax(ii) = axes('Units','normalized', ...
                'Position',[px py axw axh], ...
                'XTickLabel','', ...
                'YTickLabel',''); hold on;
            px = px+axw+gap(2);
            set( ax(ii),'ylim',yrange,'xlim',xrange,'fontsize',8,'ycolor','k','xcolor','k','color','k' );
            
            % add the y-axis to the left column
            if mod( ii,cols ) == 1
                ylabel( '\muV' );
                set( gca,'ycolor','w' );
                rowCounter = rowCounter + 1;
            end
        end
        py = py-axh-gap(1);
        %set( ax(i),'LooseInSet',get( ax(i),'TightInSet' ),'position',[.12 .08 .86 .88] );
    end
    linkaxes( get( gcf,'Children'),'xy' );
end