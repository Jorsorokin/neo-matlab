function p = st_check_inputs( inputs )
    % function p = st_check_inputs( inputs )
    %
    % parses the optional name-value inputs passed into the "sortTool" GUI
    % upon creation of the GUI instance
    names = {'data','labels','times','trials','projection'};
    defaults = nan( 1,numel( names ) );
    p = inputParser;
    for j = 1:numel( names )
        p.addParameter( names{j},defaults(j) );
    end

    % parse the inputs
    p.parse(inputs{:});
    p = p.Results;
    if max( any( p.data ) ) || max( any( p.projection ) )
        nsp = max( size( p.data,2 ),size( p.projection,1 ) );
        if isnan( p.labels )
            p.labels = zeros( 1,nsp,'uint8' );
        end
        if ~max( any( p.trials ) )
            p.trials = ones( 1,nsp,'uint8' );
        end
    end
    
    % now make into the right format
    if ~isrow( p.labels )
        p.labels = p.labels';
    end
    if ~isrow( p.times )
        p.times = p.times';
    end
    if ~isrow( p.trials )
        p.trials = p.trials';
    end
end