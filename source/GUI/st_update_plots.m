function st_update_plots( handles ) 
    % function st_update_plots( handles )
    %
    % updates all plots associated with the sortTool GUI to reflect newly
    % selected data and/or new labels. Only updates plots if they have been
    % created, AND if they are visible (no need to spend memory plotting
    % invisible lines).
    
    st_clear_plots( handles );
     
    % update the main projection plot
    if ~isempty( handles.mapplot.Children ) % i.e. if any projection is there
        handles.mapplot.Children.CData = handles.plotcolor;
    end

    % raw waveform plot (if visible / initiated)\
    if ~isa( handles.waveformplot,'double' )   
        if handles.waveformplot.isvalid && strcmp( handles.waveformplot.Visible,'on' )
            st_plot_waveforms( handles );
        end
    end

    % isi plot
    if ~isa( handles.isiplot,'double' )
        if handles.isiplot.isvalid && strcmp( handles.isiplot.Visible,'on' )
            st_plot_isi( handles );
        end
    end

    % raster plot
    if ~isa( handles.rasterplot,'double' )
        if handles.rasterplot.isvalid && strcmp( handles.rasterplot.Visible,'on' )
            st_plot_raster( handles );
        end
    end
    
    % correlogram plot
    if ~isa( handles.xcorrplot,'double' ) 
        if handles.xcorrplot.isvalid && strcmp( handles.xcorrplot.Visible,'on' )
            st_plot_correlograms( handles );
        end
    end
    
    % drift plot
    if ~isa( handles.driftplot,'double' )
        if handles.driftplot.isvalid && strcmp( handles.driftplot.Visible,'on' )
            st_plot_drift( handles );
        end
    end
    
    % dvdt plot
    if ~isa( handles.dvdtplot,'double' )
        if handles.dvdtplot.isvalid && strcmp( handles.dvdtplot.Visible,'on' )
            st_plot_dvdt( handles );
        end
    end

    % quality plot
    if ~isa( handles.qualityplot,'double' ) && handles.qualityplot.isvalid
        if all( ~isnan( handles.R.clusterQuality ) ) && any( handles.selectedPoints )
            clusts = handles.labels(handles.selectedPoints);
            clusts = ismember( handles.qualityplot.Children(2).Children.XData,clusts );
            handles.qualityplot.Children(2).Children.SizeData(clusts) = 300;
        end
    end
end