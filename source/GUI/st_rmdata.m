function handles = st_rmdata( handles,pts )
    % function handles = st_rmdata( handles,pts )
    %
    % eliminates data selected data points (in handles.selectedPoints)
    
    handles.selectedPoints(pts) = [];
    handles.projection(pts,:) = [];
    handles.labels(pts) = [];
    handles.data(:,pts,:) = [];
    handles.plotcolor(pts,:) = [];
    handles.R.keptPts(pts) = false;
    
    if ~any( isnan( handles.trials ) )
        handles.trials(pts) = [];
    end
    if ~any( isnan( handles.times ) )
        handles.times(pts) = [];
    end
    if isfield( handles.R,'probabilities' ) & ~any( isnan( handles.R.probabilities ) )
        handles.R.probabilities(pts,:) = [];
    end
    if isfield( handles,'probabilities' ) & ~any( isnan( handles.probabilities ) )
        handles.probabilities(pts,:) = [];
    end
    if ~any( isnan( handles.mask ) )
        handles.mask(:,pts) = [];
    end
end
