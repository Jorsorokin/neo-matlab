function handles = st_rightclick_delete( hObject,handles,callback );
% handles = st_rightclick_delete( hObject,handles,callback );
%
% removes a label or data point nearest to mouse rightclick from the projection plot

% get closest point to mouse
mouse = get( handles.mapplot,'CurrentPoint' );
xy = mouse(1:1:2);
ax = axis;
dx = ax(2) - ax(1);
dy = ax(4) - ax(3);

[~,closestPt] = min( sqrt( ((handles.mapplot.Children.XData - xy(1))/dx).^2...
                            + ((handles.mapplot.Children.YData - xy(2))/dx).^2 ) );

% draw a white circle around the point
pos = [handles.mapplot.Children.XData(closestPt),handles.mapplot.Children.YData(closestPt)];
circle = plot( pos(1),pos(2),'wo' );

% switch for the delete choice
switch hObject.Label
    case 'label'
        handles.labels(closestPt) = 0;
        handles.mapplot.Children.CData(closestPt,:) = handles.allPlotColor(1,:);
    case 'point'
        handles.labels(closestPt) = [];
        handles.mapplot.Children.XData(closestPt) = [];
        handles.mapplot.Children.YData(closestPt) = [];
        handles.mapplot.Children.CData(closestPt,:) = [];
        handles = st_rmdata( handles,closestPt );
end

% update the plots
delete( circle );
st_update_plots( handles );

end
