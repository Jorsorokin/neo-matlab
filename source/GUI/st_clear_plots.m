function st_clear_plots( handles )
    % function st_clear_plots( handles )
    %
    % clears all plots (waveform, isi, raster) except for the main plot in
    % the sortTool GUI
    
    % check for validity of plot handles & clear each
    if ~isa( handles.waveformplot,'double' ) && isvalid( handles.waveformplot )
        arrayfun( @cla,handles.waveformplot.Children );
    end
    
    if ~isa( handles.isiplot,'double' ) && isvalid( handles.isiplot )
        cla( handles.isiplot.Children );
    end

    if ~isa( handles.rasterplot,'double' ) && isvalid( handles.rasterplot )
        arrayfun( @cla,handles.rasterplot.Children );
    end
    
    if ~isa( handles.xcorrplot,'double' ) && isvalid( handles.xcorrplot ) 
        delete( handles.xcorrplot.Children );
    end
    
    if ~isa( handles.driftplot,'double' ) && isvalid( handles.driftplot )
        arrayfun( @cla,handles.driftplot.Children );
    end
    
    if ~isa( handles.dvdtplot,'double' ) && isvalid( handles.dvdtplot )
        cla( handles.dvdtplot.Children );
    end

    if ~isa( handles.qualityplot,'double' ) && isvalid( handles.qualityplot )
        ax = handles.qualityplot.Children(2);
        if ~isempty( ax.Children )      
            ax.Children.SizeData(:) = 50;
        end
    end
end