function handles = st_rmdata( handles,pts )
    % function handles = st_rmdata( handles,pts )
    %
    % eliminates data selected data points (in handles.selectedPoints), and
    % turns handles.model.keptPts = false for these points;
    handles.selectedPoints(pts) = [];
    handles.projection(pts,:) = [];
    handles.labels(pts) = [];
    handles.data(:,pts) = [];
    handles.model.keptPts(pts) = [];
    if ~any( isnan( handles.trials ) )
        handles.trials(pts) = [];
    end
    if ~any( isnan( handles.times ) )
        handles.times(pts) = [];
    end
    if ~isnan( handles.model.probabilities )
        handles.model.probabilities(pts,:) = [];
    end
    
end
