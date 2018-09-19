function block = add_changroup_to_block( block,chanNum,varargin )
    % block = add_changroup_to_block( block,chanNums,(distances) )
    %
    % creates a channelindx object and assigns electrodes with IDs
    % contained in "chanNum", assuming electrodes have already been added
    % to the block
    
    % pull out the appropriate electrodes
    electrodeID = [block.getChild( 'Electrode' ).electrodeNum];
    electrode = block.getChild( 'Electrode',find( ismember( electrodeID,chanNum ) ) );
    if isempty( electrode )
        fprintf( 'No electrodes for this group were found in the block\n' );
        return
    end
    
    % create our ChannelIndex object
    channelindex = ChannelIndex( block.nChanInds + 1 );
    channelindex.channels = chanNum; % i.e. relative to list of all electrodes in block parent 
    channelindex.addChild( electrode );
    
    % add the electrodes and distances if provided
    if nargin > 2 && ~isempty( varargin{1} )
        channelindex.chanDistances = varargin{1};        
        [~,channelindex.chanMap] = sort( channelindex.chanDistances,'descend' );
    
        if size( channelindex.chanMap,2 ) > 1
            channelindex.chanMap(:,1) = [];
        end
    end
    
    block.addChild( channelindex );
end