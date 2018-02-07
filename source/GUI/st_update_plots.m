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
        children.CData = handles.plotcolor;
    end

    % raw waveform plot (if visible / initiated)\
    if ~isa( handles.waveformplot(1),'double' )   
        if handles.waveformplot(1).isvalid && strcmp( handles.waveformplot(1).Parent.Visible,'on' )
            st_plot_waveforms( handles );
        end
    end

    % isi plot
    if ~isa( handles.isiplot,'double' )
        if handles.isiplot.isvalid && strcmp( handles.isiplot.Parent.Visible,'on' )
            st_plot_isi( handles );
        end
    end

    % raster plot
    if ~isa( handles.rasterplot,'double' )
        if handles.rasterplot.isvalid && strcmp( handles.rasterplot.Parent.Visible,'on' )
            st_plot_raster( handles );
        end
    end

    % quality plot
    if ~isa( handles.qualityplot,'double' ) && handles.qualityplot.isvalid
        if all( ~isnan( handles.R.clusterQuality ) ) && any( handles.selectedPoints )
            clusts = handles.labels(handles.selectedPoints);
            ax = handles.qualityplot.Children(2);
            clusts = ismember( ax.Children.XData,clusts );
            ax.Children.SizeData(clusts) = 300;
        end
    end
end