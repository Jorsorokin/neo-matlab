function st_update_plots( handles ) 
    % function st_update_plots( handles )
    %
    % updates all plots associated with the sortTool GUI to reflect newly
    % selected data and/or new labels. Only updates plots if they have been
    % created, AND if they are visible (no need to spend memory plotting
    % invisible lines).
    
    % first, clear the waveform, isi, and raster plots
    st_clear_plots( handles );
     
    % update the main projection plot
    children = handles.mapplot.Children;
    if ~isempty( children ) % i.e. if any projection is there
        if numel( children ) > 1 && ~isa( handles.selectedPointPlot,'double' )
            mainPlot = ~ismember( children,handles.selectedPointPlot );
        else
            mainPlot = 1;
        end
        children(mainPlot).CData = handles.plotcolor;
    end

    % raw waveform plot (if visible / initiated)\
    if ~isa( handles.waveformplot(1),'double' )   
        if isvalid( handles.waveformplot(1) ) && strcmp( handles.waveformplot(1).Parent.Visible,'on' )
            st_plot_waveforms( handles );
        end
    end

    % isi plot
    if ~isa( handles.isiplot,'double' )
        if isvalid( handles.isiplot ) && strcmp( handles.isiplot.Parent.Visible,'on' )
            st_plot_isi( handles );
        end
    end

    % raster plot
    if ~isa( handles.rasterplot,'double' )
        if isvalid( handles.rasterplot ) && strcmp( handles.rasterplot.Parent.Visible,'on' )
            st_plot_raster( handles );
        end
    end
end