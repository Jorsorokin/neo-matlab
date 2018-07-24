function sort_blocks( files,templates,chindex,thresh )
    % sort_blocks( files,templates,chindex,thresh )
    %
    % sorts the channelindex specified by 'chindex' in for each block file
    % contained in the cell array "files" using the templates structure
    % "templates". Doesn't resolve overlap and assumes masks / channel maps
    % are used in the template matching step
    
    outdir = dir( files{1} );
    outdir = outdir.folder; 
    
    for f = 1:numel( files )
        load( files{f} );
        chind = block.getChild( 'ChannelIndex',chindex );
        chind.mergeNeurons( [chind.getChild( 'Neuron' ).ID] );
        chind.sortSpikes_using_templates( templates,true,true,false,thresh );
        block.write( outdir );
    end
end