function ax = st_create_waveplot( handles )
    % ax = st_create_waveplot( handles )
    % creates an instance of the "waveformplot" used for the "sortTool" GUI
    %
    % handles is a structure with variables referenced by the main GUI
    figure( 'Visible','off','color','k','position',[4,530,340,250] );
    set(gcf,'Name','Raw waveforms');
    if ~isnan( handles.data )
        yrange = [min( min( handles.data ) ),max( max( handles.data ) )];
        xrange = [1, size( handles.data,1 )];
    else
        yrange = [-500 200];
        xrange = [0 1];
    end
    ax = subplot( 1,1,1 ); hold on;
    set( ax,'tickdir','out','box','off','color','k',...
        'ylim',yrange,'xlim',xrange,'fontsize',8,'xcolor','k','ycolor','w' );
    set( ax,'LooseInSet',get( ax,'TightInSet' ),'position',[.12 .08 .86 .88] ); 
    ylabel('\muV');
end