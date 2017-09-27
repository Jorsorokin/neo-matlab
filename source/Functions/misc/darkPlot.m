function darkPlot( fig )
% darkPlot( fig )
%
% reformats the figure specified by the figure handle 
% "fig" to a dark-theme

set( fig,'color','k' );
ax = get( fig,'Children' );
for i = 1:numel( ax )
    set( ax(i),'color','k','xcolor','w','ycolor','w',...
        'tickdir','out','box','off' );
    set( get( ax(i),'Title' ),'color','w' );
end

end

