function plot_spikes_over_data( data,fs,sptm,mask,varargin )
    % plot_spikes_over_data( data,fs,sptm,mask,(labels,linesep) )
    %
    % plots markers for spikes overlaid on raw data. Can optionaly specify labels to
    % plot the markers with different colors
    
    [~,m] = size( data );
    
    if nargin > 4 
        labels = varargin{1};
    else
        labels = ones( 1,numel( sptm ),'int8' );
    end
    
    if nargin > 5 
        linesep = varargin{2};
    else
        linesep = [];
    end
    
    [fh,lh] = multisignalplot( data,fs,'w',linesep ); hold on;
    dY = abs( mean( lh(1).YData ) - mean( lh(2).YData ) ) / 4;
    fh.Visible = 'off';
    
    for ch = 1:m
        idx = mask(ch,:) > 0;
        s = scatter( sptm(idx),ones( 1,nnz( idx ) )*mean( lh(ch).YData )+dY,[],labels(idx) );
        set( s,'Marker','x','SizeData',40,'LineWidth',2 );
        %lh(ch).Color = 'w';
    end
    
    darkPlot( gcf );
    fh.Visible = 'on';
end
    
    
    
    
    