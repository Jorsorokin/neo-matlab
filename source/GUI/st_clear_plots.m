function st_clear_plots( handles )
    % function clear_plots( handles )
    %
    % clears all plots (waveform, isi, raster) except for the main plot in
    % the sortTool GUI
    
    % check for validity of plot handles & clear each
    if ~isa( handles.waveformplot,'double' ) && isvalid( handles.waveformplot )
        cla( handles.waveformplot );
    end
    if ~isa( handles.isiplot,'double' ) && isvalid( handles.isiplot )
        cla( handles.isiplot );
    end
    if ~isa( handles.rasterplot,'double' ) && isvalid( handles.rasterplot )
        cla( handles.rasterplot );
    end
end