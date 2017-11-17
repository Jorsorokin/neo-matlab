function st_clear_plots( handles )
    % function st_clear_plots( handles )
    %
    % clears all plots (waveform, isi, raster) except for the main plot in
    % the sortTool GUI
    
    % check for validity of plot handles & clear each
    if ~isa( handles.waveformplot(1),'double' ) && isvalid( handles.waveformplot(1) )
        for c = 1:numel( handles.waveformplot )
            cla( handles.waveformplot(c) );
        end
    end
    
    if ~isa( handles.isiplot,'double' ) && isvalid( handles.isiplot )
        cla( handles.isiplot );
    end

    if ~isa( handles.rasterplot,'double' ) && isvalid( handles.rasterplot )
        cla( handles.rasterplot );
    end

    if ~isa( handles.qualityplot,'double' ) && isvalid( handles.qualityplot )
        ax = handles.qualityplot.Children(2);
        if ~isempty( ax.Children )      
            ax.Children.SizeData(:) = 50;
        end
    end
end