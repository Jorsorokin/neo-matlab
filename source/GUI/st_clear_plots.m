function st_clear_plots( handles )
    % function clear_plots( handles )
    %
    % clears all plots (waveform, isi, raster) except for the main plot in
    % the sortTool GUI
    
    % check for validity of plot handles & clear each
    for c = 1:numel( handles.waveformplot )
        if ~isa( handles.waveformplot(c),'double' ) && isvalid( handles.waveformplot(c) )
            cla( handles.waveformplot(c) );
        end
    end
    if ~isa( handles.isiplot,'double' ) && isvalid( handles.isiplot )
        cla( handles.isiplot );
    end
    if ~isa( handles.rasterplot,'double' ) && isvalid( handles.rasterplot )
        cla( handles.rasterplot );
    end
end