function varargout = sortTool(varargin)
    % code for sortTool.fig
    %
    % GUI for facilitating spike sorting of large databases. Inputs are optional,
    % as one can load in the various files through the GUI itself. This GUI
    % works equally well for any data that needs sorting (genomic, EPSPs,
    % network graphs, etc.). The only difference is that some features may
    % be irrelevant (i.e. ISI and rasters, for instance).
    %
    % The GUI projects high-dimensional data using a variety of projection
    % methods that the user can specify. These projection methods are a
    % mixture of those written by the "FastICA" team, Laurens van der
    % Maaten's wonderful dimension reduction toolbox, and some written by
    % myself. The user can choose how many dimensions to project onto, 
    % as well as which two dimensions (i.e. which columns of the projected
    % data matrix) to visualize. Clicking "re-plot" will plot the
    % projections using a new viewpoint if specified.
    %
    % Following projections, one can choose to automatically sort
    % the clusters - and can further specify the number of clusters to sort
    % or allow an algorithm to solve for this number - using either
    % Expectation Maximization for Gaussian Mixture Models (EM-GMM),
    % Expectation Maximization for T-distribution Mixture Models (EM-TMM),
    % K-means (Km), or Variational Bayes (VB). 
    %
    % Finally, one can also manually delete, merge, and create new clusters 
    % by clicking on the "Manual sort" button, and can select individual 
    % data points or individual clusters for viewing associated waveforms, 
    % ISI, and raster plots by clicking the "Select points" and the
    % "Select cluster" buttons, respectively. To remove all cluster labels, 
    % click "Delete labels", and to close the GUI and extract the final
    % cluster labels, click "Finish". 
    %
    % INPUTS:
    %   ( data ) - the raw waveforms in column-order (each column is one
    %              observation, rows are variables). Can be 2- or 3-dimensional
    %
    %   ( times ) - the times of the waveforms as one long vector, in seconds
    %
    %   ( projection ) - the projections of the data onto lower dims
    %
    %   ( labels ) - the labels from a previous sorting routine
    %
    %   ( trials ) - trial-labels of the waveforms (i.e. one can upload a
    %                a large data matrix with waveforms from different trials,
    %                but to get an accurate representation of ISI one needs
    %                information regarding which spikes are from which
    %                trial. Default is to assume the same trail).
    %
    %   ( mask ) - a sparse or full c x m matrix, with 0 <= (i,j) <= 1. The
    %              mask matrix defines which channels 1:c are "involved"
    %              with spikes 1:m. A mask matrix is used when computing 
    %              features by "masking" channels not involved with any spike.
    %
    %   ( location ) - a c x q matrix specifying the locations of each electrode.
    %                  Useful for dense multi-electrode probe configurations. Note,
    %                  this is only used if "data" is 3-dimensional (with 3rd dimension = channels)
    % 
    % OUTPUTS:
    %   ( labels ) - the final labels assigned to each spike. Note, a label
    %                of "0" represents unsorted spikes.
    %
    %   ( projection ) - the projections of the spike waveforms onto a
    %                    lower-dimensional subspace specified in the GUI
    %
    %   ( model ) - a structure containing information on the projection /
    %               sorting including 
    %                   mapping         - the projection mapping structure (if applicable)
    %                   projMethod      - the projection method (PCA, ICA, LE,...)
    %                   sortMethod      - the sorting method (Km, EM-GMM, DBSCAN,... )
    %                   sortModel       - the sorting model structure (if applicable)
    %                   probabilities   - the probabilities of each spike
    %                                       belonging to each kth cluster
    %
    % Last edited by Jordan Sorokin, 9/8/2017

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @sortTool_OpeningFcn, ...%'gui_OutputFcn',  @sortTool_OutputFcn,
                       'gui_OutputFcn',  @sortTool_OutputFcn,...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    % End initialization code - DO NOT EDIT

    
    % --- Executes just before sortTool is made visible.
    function sortTool_OpeningFcn(hObject, eventdata, handles, varargin)
        
        % Choose default command line output for sortTool
        handles.output = hObject;

        % set up the figure
        set( 0,'defaultfigurewindowstyle','normal' );
        handles.figure1.Name = 'Dimensionality Reduction & Sorting Tool';
        handles.figure1.Visible = 'off';                % turns of the fig temporarily
        set( handles.figure1,'toolbar','figure' );      % turns on the toolbar
        figureMenu = uimenu( 'Text','Load' );           % creates a 'load' menu button
        handles.figureMenu = uimenu( figureMenu,...     % creats variables to load
            'Text',{'data','times','projection',...
            'labels','trials','mask','location'}...
            'Callback',@loaddata_Callback );
        set( handles.figure1,'WindowButtonMotionFcn',@mouseXY); % allows for tracking location of mouse

        % set up user-defined variables in the handles structure
        handles.data = nan;                             % the waveforms
        handles.availableData = false;                  % to avoid "gathering" gpuArray for each logical check
        handles.times = nan;                            % for plotting rasters/ISI
        handles.projection = nan;                       % the projection using the "projectMethod"    
        handles.availableProjection = false;            % to avoid "gathering" gpuArray for each logical check
        handles.labels = nan;                           % previous sorting labels
        handles.trials = nan;                           % n-dim vector used for plotting rasters
        handles.mask = nan;                             % mask vector (see "double_flood_fill.m")
        handles.location = nan;                         % vector specifying (x,y) location of each electrode
        handles.selectedPoints = nan;                   % points highlighted on the projection plot
        handles.manualClust = {};                       % will be added to while with manual cluster assignment
        handles.plotDims = [1,2];                       % default plotting 1st & 2nd columns
        handles.nDim = 3;                               % defaults to 3D projection 
        
        % create our model structure for ouput
        handles.model = struct; 
        handles.model.projectMethod = 'PCA';            % defaults to PCA first
        handles.model.sortMethod = 'EM-GMM';            % defaults to EM gaussian mixture modeling
        handles.model.mapping = nan;                    % the projection mapping
        handles.model.keptPts = nan;                    % vector indicating which points we manually delete/remove labels

        % plotting colors for different labels
        handles.allPlotColor = [ [.75, .75, .75];...    % light gray        (0)
                                 [0.1, 0.74, 0.95];...  % deep sky-blue     (1)
                                 [0.95, 0.88, 0.05];... % gold/yellow       (2)
                                 [0.80, 0.05, 0.78];... % magenta           (3)
                                 [0.3, 0.8, 0.20];...   % lime green        (4)
                                 [0.95, 0.1, 0.1];...   % crimson red       (5)   
                                 [0.64, 0.18, 0.93];... % blue-violet       (6)
                                 [0.88, 0.56, 0];...    % orange            (7)
                                 [0.4, 1.0, 0.7];...    % aquamarine        (8)
                                 [0.95, 0.88, 0.7];...  % salmon-yellow     (9)
                                 [0, 0.2, 1];...        % blue              (10)
                                 [1, 0.41, 0.7];...     % hot pink          (11)
                                 [0.5, 1, 0];...        % chartreuse        (12)
                                 [0.6, 0.39, 0.8];...   % amtheyist         (13)
                                 [0.82, 0.36, 0.36,];...% indian red        (14)
                                 [0.53, 0.8, 0.98];...  % light sky blue    (15)
                                 [0, 0.6, 0.1];...      % forest green      (16)
                                 [0.65, 0.95, 0.5];...  % light green       (17)
                                 [0.85, 0.6, 0.88];...  % light purple      (18)
                                 [0.90, 0.7, 0.7];...   % light red         (19)
                                 [0.1, 0.1, 0.5];...    % dark blue         (20)
                                ];
                                
        % repeat colors in case we have many identified neurons
        handles.allPlotColor = [handles.allPlotColor; repmat( handles.allPlotColor(2:end,:),4,1 )]; 
        handles.plotcolor = nan;    
        handles.buttonPressColor = [0.5, 0.8, 0.95];
        
        % GET THE OPTIONAL INPUTS
        p = st_check_inputs( varargin );
        handles.data = p.data;
        if any( ~isnan( gather( handles.data ) ) )
            handles.availableData = true;
        end
        handles.projection = p.projection;
        if any( ~isnan( gather( handles.projection ) ) )
            handles.availableProjection = true;
        end
        handles.labels = p.labels;
        handles.times = p.times;
        handles.trials = p.trials;
        handles.mask = p.mask;
        handles.location = p.location;
        clear p
               
        % get the selectedPoints if data available
        if ~isnan( handles.labels )
            handles.selectedPoints = false( 1,numel( handles.labels ) );
            handles.model.keptPts = single( 1:numel( handles.labels ) );
        end
        
        % INITIATE THE PLOTS
        axes( handles.mapplot ); hold on;
        set( gca,'tickdir','out','box','off','fontsize',8,...
            'xcolor','k','ycolor','k','ylim',[-1 1],'xlim',[-1 1]);
        handles.waveformplot = nan;
        handles.isiplot = nan;
        handles.rasterplot = nan;
        handles.plotcolor = update_plot_colors( handles );
        
        % plot the projections if provided
        if handles.availableProjection
            st_plot_projections( handles );
        end
        
        % turn on the main figure
        handles.figure1.Visible = 'on';
        
        % Update handles structure
        guidata( hObject, handles );
        
        % keep track of the GUI visibility before running the OutputFcn
        waitfor( handles.figure1,'Visible','off' );
        

    % --- Outputs from this function are returned to the command line.
    function varargout = sortTool_OutputFcn(hObject, eventdata, handles)
        % varargout  cell array for returning output args (see VARARGOUT);
        varargout{1} = handles.labels;
        varargout{2} = handles.projection;
        varargout{3} = handles.model;
        
        close( handles.figure1 );
        if ~isa( handles.waveformplot(1),'double' ) && isvalid( handles.waveformplot(1) )
            close( handles.waveformplot(1).Parent );
        end
        if ~isa( handles.isiplot,'double' ) && isvalid( handles.isiplot )
            close( handles.isiplot.Parent );
        end
        if ~isa( handles.rasterplot,'double' ) && isvalid( handles.rasterplot )
            close( handles.rasterplot.Parent );
        end
        
                 
    % --- clears the current figure
    function CloseMenuItem_Callback(hObject, eventdata, handles)
        selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                             ['Close ' get(handles.figure1,'Name') '...'],...
                             'Yes','No','Yes');
        if strcmp(selection,'No')
            return;
        end
        delete( handles.figure1 ); 
        close all
        

    % --- closes GUI and returns outputs
    function savebutton_Callback(hObject, eventdata, handles)
        % closes the GUI   
        handles.figure1.Visible = 'off';
        
        
    % --- Deletes labels assigned to points
    function deletebutton_Callback(hObject,eventdata,handles)
        if isnan( handles.labels )
            return;
        end
        
        % clear previous labels, update the plots
        handles = st_clearSelectedData( handles );
        handles.labels(:) = 0;
        handles.plotcolor = update_plot_colors( handles );
        st_update_plots( handles );
        guidata( hObject,handles );


    % --- gets current position of mouse and highlights pushbuttons if necessary
    function mouseXY( gcbo,eventdata,handles )
        mousePosition = get( handles.figure1, 'CurrentPoint' );
        if any( mousePosition ~= handles.mousePosition ) % we've moved the mouse
            handles.mousePosition = mousePosition;
        else
            return % no movement, so don't loop over the pushbuttons
        end

        % search over all pushbutton objects
        buttons = findjobj( handles.figure1,'Type','pushbutton' );
        for j = 1:numel( buttons )
            pos = buttons(j).getpixelposition( buttons(j),'true' ); % relative to figure
            inX = any( ismember( round( pos(1):pos(1)+pos(3) ),round( handles.mousePosition(1) ) ) );
            inY = any( ismember( round( pos(2):pos(2)+pos(4) ),round( handles.mousePosition(2) ) ) );
            if inX & inY
                buttons(j).BackgroundColor = [0.4 0.6 0.8]; % light blue
            end
        end

        guidata( handles.figure1, handles );


% loading & projections
% =======================================================
    % --- Executes on button press in menu item 'load'. Allows user to load in data into GUI
    function loaddata_Callback(hObject, eventdata, handles)

        st_clear_plots( handles )

        % load in the variable names by creating a table of variable names
        vars = who;
        [ind,val] = listdlg( 'PromptString','Select variable',...
                        'SelectionMode','Single',...
                        'ListString',vars );
        if val == 0
            return
        else
            % get the name of the variable
            name = vars{ind};

            % load the variable into the handles object
            x = eval( 'base',name ); % loads in data into 'x'
            switch hObject.Text
                case 'data'
                    handles.data = x;
                    handles.availableData = true;
                    handles.labels = zeros( 1,size( x,2 ),'uint8' );
                case 'times'
                    handles.times = x;
                case 'projection'
                    handles.projection = x;
                    handles.availableProjection = true;
                case 'labels'
                    handles.labels = x;
                case 'trials'
                    handles.trials = x;
                case 'mask'
                    handles.mask = x;
                case 'location'
                    handles.location = x;
            end
        end

        % update the GUI
        guidata( hObject, handles );
        
        
    % --- Executes during object creation, after setting all properties.
    function mappingmethod_CreateFcn(hObject, eventdata, handles)
        % establishes the possible mapping methods in the dropdown bar

        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
             set(hObject,'BackgroundColor','white');
        end

        % set the available options
        set( hObject, 'String',...
            {'PCA (principal components analysis)',...
            'kPCA (kernel PCA)',...
            'pPCA (probabilistic PCA',...
            'NCA (neighborhood components analysis)',...
            'LPP (locality preserving projection)',...
            'ICA (independent components analysis',...
            'LLE (local linear embedding)',...
            'HLLE (hessian local linear embedding)',...
            'NPE (neighborhood preserving embedding)',...
            'SPE (stochastic proximity embedding)',...
            'MVU (maximum variance unfolding)',...
            'FastMVU',...
            'LE (laplacian eigenmaps)',...
            'DM (diffusion maps)',...
            'SNE (stochastic neighborhood embedding)',...
            'tSNE (t-dist. stochastic neighborhood embedding)',...
            'Isomap',...
            'Sammon',...
            'Autoencoder'} );

        % update the GUI
        guidata( hObject, handles );   


    % --- Executes on selection change in mappingmethod.
    function mappingmethod_Callback(hObject, eventdata, handles)
        % selects a new mapping method

        allmethods = {'PCA','KernelPCA','ProbPCA','NCA','LPP','ICA','LLE',...
                        'HessianLLE','NPE','SPE','MVU','FastMVU','Laplacian',...
                        'DiffusionMaps','SNE','tSNE','Isomap','Sammon','Autoencoder'};
        value = get( hObject,'value' );
        handles.model.projectMethod = allmethods{value};

        % update the GUI
        guidata( hObject, handles );

        
    % --- Executes on button press in projectdata.
    function projectdata_Callback(hObject,eventdata,handles)
        % projects waveforms using the chosen mapping method

        % do the projections
        [handles.projection,handles.model.mapping] = st_project_data( handles );
        handles.availableProjection = true;
        
        % now plot onto the main figure
        handles.plotDims( handles.plotDims > handles.nDim ) = handles.nDim;
        handles.mapplotMenu = st_plot_projections( handles );
        handles.selectedPointPlot = st_plotSelectedData( handles );

        % update the GUI
        guidata( hObject, handles );

                
    % --- Executes during object creation, after setting all properties.
    function dim1_CreateFcn(hObject, eventdata, handles) 

        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end
        set(hObject,'String',1); % automatically sets to 1 when initializing
        guidata( hObject,handles );

        
    function dim1_Callback(hObject, eventdata, handles)
        % specifies the first dimension to display

        handles.plotDims(1) = str2double( get( hObject,'String' ) );
        guidata( hObject,handles );

        
    % --- Executes during object creation, after setting all properties.
    function dim2_CreateFcn(hObject, eventdata, handles)

        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end

        set(hObject,'String',2); % automatically sets to 2 when initializing
        guidata( hObject, handles );
        
        
    function dim2_Callback(hObject, eventdata, handles)
        % specifies the second dimension to display
        
        handles.plotDims(2) = str2double( get( hObject,'String' ) );
        guidata( hObject,handles );     
        

    % --- Executes during object creation, after setting all properties.          
    function nDim_CreateFcn(hObject, eventdata, handles)

        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end
        set( hObject,'String',3 ); % default = 3
        guidata( hObject,handles );
        

    function nDim_Callback(hObject, eventdata, handles)
        % specifies the number of dimensions to project onto using the specified mapping method
  
        handles.nDim = str2double( get( hObject,'String' ) );
        guidata( hObject,handles );    
        

% plotting procedures (updating colors, replotting, etc.)
% ======================================================= 
     % --- Executes on button press in replogt
    function replogt_Callback(hObject, eventdata, handles)
        % update the plots according to new labels or projection dims

        handles.mapplotMenu = st_plot_projections( handles );
        handles.selectedPointPlot = st_plotSelectedData( handles );
        st_update_plots( handles );
        guidata( hObject,handles );
        
        
    % --- Executes on button press in selectdata
    function selectdata_Callback(hObject, eventdata, handles)
        % allows user to select individual data points in the main plot  

        if ~handles.availableData
            return
        end

        % get the selection value, plot the selected points
        pts = st_selectDataLasso( handles.mapplot );
        if iscell( pts )
            pts = pts{end};
        end
        handles.selectedPoints(pts) = true;
        handles.selectedPointPlot = st_plotSelectedData( handles );
        
        % update 
        st_update_plots( handles );
        guidata( hObject,handles );
        
    
    % --- executes on select cluster button press
    function selectclust_Callback(hObject,eventdata,handles)
        % allows user to select individual clusters in the main plot

        if ~handles.availableData
            return
        end
        
        % get the selected point and those belonging to same cluster
        pts = st_selectDataCluster( handles.mapplot,handles.labels );
        if iscell( pts )
            pts = pts{end};
        end
        handles.selectedPoints(pts) = true;
        handles.selectedPointPlot = st_plotSelectedData( handles );
        
        % update
        st_update_plots( handles );
        guidata( hObject,handles );
        

    % --- clears the manually-selected data on button press
    function cleardata_Callback(hObject, eventdata, handles)

        % clear all data selection
        handles = st_clearSelectedData( handles );
        
        % update
        guidata( hObject,handles );
        
                    
    % --- updates plotting colors of all points
    function colors = update_plot_colors( handles )
        
        % check if any labels (i.e. any data)
        if any( isnan( handles.labels ) )
            colors = nan;
            return
        else
            colors = zeros( numel( handles.labels ),3,'single' );
        end
        
        % update plotting colors for each data point
        uID = unique( handles.labels );
        for i = uID
            idx = handles.labels==i;
            colors(idx,:) = repmat( single( handles.allPlotColor(i+1,:) ),sum( idx ),1 ); % since uID==0 is first color
        end
            
               
    % --- executes with plotspikes push (allows waveform plotting)
    function plotspikes_Callback(hObject, eventdata, handles)
        
        % create a new figure if needed
        if isa( handles.waveformplot(1),'double' ) || ~isvalid( handles.waveformplot(1) )
            handles.waveformplot = st_create_waveplot( handles ); 
        end
            
        % plot
        switch handles.waveformplot(1).Parent.Visible
            case 'on'
                cla( handles.waveformplot );
                handles.waveformplot(1).Parent.Visible = 'off'; % hide the plot
            otherwise
                handles.waveformplot(1).Parent.Visible = 'on';
                st_plot_waveforms( handles )
        end
        
        % update
        guidata( hObject,handles );    
        
        
     % --- executes with plotspikes push (allows waveform plotting)
    function plotisi_Callback(hObject, eventdata, handles)
        
        % create a new figure if needed
        if isa( handles.isiplot,'double') || ~isvalid( handles.isiplot )
            handles.isiplot = st_create_isiplot( handles );
        end
        
        % check if visible; if not, plot
        switch handles.isiplot.Parent.Visible
            case 'on'
                cla( handles.isiplot );
                handles.isiplot.Parent.Visible = 'off';
            otherwise
                st_plot_isi( handles );
        end
        
        % update
        guidata( hObject,handles );  
        
        
    % --- executes with plotspikes push (allows waveform plotting)
    function plotraster_Callback(hObject, eventdata, handles)
        
        % create a new figure if needed
        if isa( handles.rasterplot,'double' ) || ~isvalid( handles.rasterplot )
            handles.rasterplot = st_create_rasterplot( handles );
        end
        
        % check if Visible; if not, plot the raster
        switch handles.rasterplot.Parent.Visible
            case 'on'
                cla( handles.rasterplot );
                handles.rasterplot.Parent.Visible = 'off';
            otherwise
                st_plot_raster( handles );
        end
        
        % update
        guidata( hObject,handles );  

    % --- executes with right-click 'assign' (re-assigns single point label)
    function rightclick_assign( hObject,handles,callback );

        % get closest point to mouse
        mouse = get( handles.mapplot,'CurrentPoint' );
        xy = mouse(1:1:2);
        ax = axis;
        dx = ax(2) - ax(1);
        dy = ax(4) - ax(3);

        [~,closestPt] = min( sqrt( ((handles.mapplot.Children.XData - xy(1))/dx).^2...
                                    + ((handles.mapplot.Children.YData - xy(2))/dx).^2 ) );

        % draw a white circle around the point
        pos = [handles.mapplot.Children.XData(closestPt),handles.mapplot.Children.YData(closestPt)];
        circle = plot( pos(1),pos(2),'wo' );

        % switch for the delete choice
        switch hObject.Label
            case 'new label'
                handles.labels(closestPt) = max( handles.labels ) + 1;
                handles.mapplot.Children.CData(closestPt,:) = handles.allPlotColor(handles.labels(closestPt,:));
            case 'existing label'
                disp( 'currently under development' );
        end

        % update the plots
        delete( circle );
        st_update_plots( handles );
        guidata( hObject,handles );


    % --- executes with right-click 'delete' (deletes single label or point)
    function rightclick_delete( hObject,handles,callback );

        % get closest point to mouse
        mouse = get( handles.mapplot,'CurrentPoint' );
        xy = mouse(1:1:2);
        ax = axis;
        dx = ax(2) - ax(1);
        dy = ax(4) - ax(3);

        [~,closestPt] = min( sqrt( ((handles.mapplot.Children.XData - xy(1))/dx).^2...
                                    + ((handles.mapplot.Children.YData - xy(2))/dx).^2 ) );

        % draw a white circle around the point
        pos = [handles.mapplot.Children.XData(closestPt),handles.mapplot.Children.YData(closestPt)];
        circle = plot( pos(1),pos(2),'wo' );

        % switch for the delete choice
        switch hObject.Label
            case 'label'
                handles.labels(closestPt) = 0;
                handles.mapplot.Children.CData(closestPt,:) = handles.allPlotColor(1,:);
            case 'point'
                handles.labels(closestPt) = [];
                handles.mapplot.Children.XData(closestPt) = [];
                handles.mapplot.Children.YData(closestPt) = [];
                handles.mapplot.Children.CData(closestPt,:) = [];
                handles = st_rmdata( handles,closestPt );
        end

        % update the plots
        delete( circle );
        st_update_plots( handles );
        guidata( hObject,handles );


% auto/manual sorting
% =======================================================
    % --- Initiates the sort-method dropdown list
    function sortmethod_CreateFcn(hObject,eventdata,handles)
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
             set(hObject,'BackgroundColor','white');
        end

        % set the available options
        set( hObject, 'String',...
            {'EM-GMM (gaussian mixture model)',...
            'mEM-GMM (masked EM-GMM)',...
            'EM-TMM (t-dist. mixture model)',...
            'Km (K-means)',...
            'VB (variational bayes)',...
            'mVB (masked VB)',...
            'DBSCAN',...
            'Spectral clustering'} );

        % update the GUI
        guidata( hObject, handles ); 
        
        
    % --- Executes on button press in sortmethod
    function sortmethod_Callback(hObject,eventdata,handles)
        % select a new sorting method

        allmethods = {'EM-GMM','mEM-GMM','EM-TMM','Km',...
            'VB','mVB','DBSCAN','Spectral'};
        value = get( hObject,'value' );
        handles.model.sortMethod = allmethods{value};

        % update the GUI
        guidata( hObject, handles );
    
        
    % --- Executes on button press in autosort.
    function autosort_Callback(hObject, eventdata, handles)        
        % automatically sort the data projections via the specified sort method

        warning off;
        
        % first project the data 
        if ~handles.availableData
            errordlg( 'No data available. Please load data!' );
            return
        elseif ~handles.availableProjection
            disp( 'Please project your data first!' );
            return
        end

        % get the selected points, if any
        if any( handles.selectedPoints )
            pts = handles.selectedPoints;
        else
            pts = ~handles.selectedPoints;
        end
        
        % request for the appropriate cluster number for parametric
        % clustering algorithms
        switch handles.model.sortMethod
            case {'mEM-GMM','EM-GMM','EM-TMM','VB','mVB','Km'}
                searchForK = questdlg( 'Search for number of clusters?' );
                switch searchForK
                    case 'Yes'
                        [~,K] = findClustNum( handles.projection(pts,:)',2,50 );
                    case 'No'
                        K = [];
                    case 'Cancel'
                        return
                end
        end
        
        % for masked EM compute the feature matrix without concatenation
        switch handles.model.sortMethod
            case {'mVB','mEM-GMM'}
                if any( isnan( handles.mask ) )
                    disp( 'must include a masking matrix for masked-EM and masked-VB sorting' );
                    return
                end

                % compute the projections
                projections = compute_spike_features( gather( permute( handles.data(:,pts,:),[2,1,3] ) ),...
                    handles.nDim,handles.model.projectMethod,handles.mask,[],false );

                % compute the expanded mask
                mask = zeros( size( projections ) );
                for c = 1:size( handles.mask,1 )
                    inds = c*handles.nDim-handles.nDim+1:c*handles.nDim;
                    mask(:,inds) = repmat( handles.mask(c,pts)',1,handles.nDim );
                end 
            otherwise
                projections = handles.projection(pts,:);
                mask = [];
        end
        
        % now do the actual clustering using the method specified
        [labels,~] = sort_clusters( projections,K,handles.model.sortMethod,mask );
        projections = gather( handles.projection(pts,:) );
        [mu,sigma] = get_cluster_description( projections,labels );
        [probs,prior] = compute_cluster_probabilities( projections,labels,mu,sigma );
        
        % change the labels to contain the original label plus consecutive
        % labels for each new cluster
        uID = unique( handles.labels );
        if any( uID > 0 )
            if ~any( handles.selectedPoints )
                labels = labels + min( uID(uID > 0) )-1; % subtract 1 to start with the minimum of the previous label
            else
                uID = unique( handles.labels(~pts) );
                uID(uID==0) = [];
                changeLabels = ismember( labels,uID );
                if ~isempty( changeLabels ) && any( changeLabels )
                    labels(changeLabels) = labels(changeLabels) - min( labels(changeLabels) ) + max( uID ) + 1; % creates new labels to avoid overlap 
                end
            end
        end
        
        % refine the clusters of this sorting run
        labels = st_refine_cluster( labels,probs,0.01 ); 
                
        % update all probabilities and labels in the handles struct
        handles.labels(pts) = labels;
        handles = update_labels( handles );
        handles.plotcolor = update_plot_colors( handles );
        st_update_plots( handles );
        guidata( hObject,handles );
        
        warning on;
    
        
    % --- Executes on button press in manualsort.
    function manualsort_Callback(hObject, eventdata, handles)
        % initiates the functionality for manual sorting (cluster cutting)

        % clear any data selection
        handles = st_clearSelectedData( handles );
        
        % check if data has been plotted
        axes( handles.mapplot );
        if isempty( handles.mapplot.Children )
            if handles.availableProjection
                handles.mapplotMenu = st_plot_projections( handles );
            end
        end        
        
        % set visibility of the new buttons to "on"
        set( handles.newLabel,'enable','on' );
        set( handles.mergeLabel,'enable','on' );
        set( handles.deleteLabel,'enable','on' );
        set( handles.finishManual,'enable','on' );
        set( handles.autosort,'enable','off' );
        set( handles.projectdata,'enable','off' );
        set( handles.loaddata,'enable','off' );
        set( handles.replogt,'enable','off' );
        
        % update
        guidata( hObject,handles );
        
        
    % --- Executes on button press in newLabel.
    function newLabel_Callback(hObject, eventdata, handles)
        % manually select a new label by lasso-ing around data points

        disp( 'Select a new cluster' );
        axes( handles.mapplot );

        % draw lasso
        clust = imfreehand( 'Closed',1 );
        nClust = sum( ~cellfun( @isempty,handles.manualClust ) );
        handles.manualClust{nClust+1} = clust;
        
        % update
        guidata( hObject,handles );

        
    % --- Executes on button press in mergeLabel.
    function mergeLabel_Callback(hObject, eventdata, handles)
        % merges two clusters together 

        disp( 'Click on two clusters to merge' );
        axes( handles.mapplot );
        
        % get the points and associated labels with each cluster
        [pts1,label1] = st_selectDataCluster( handles.mapplot,handles.labels );
        [pts2,label2] = st_selectDataCluster( handles.mapplot,handles.labels );
        handles.selectedPoints(pts1) = true;
        handles.selectedPoints(pts2) = true;
        handles.selectedPointPlot = st_plotSelectedData( handles );
        merge = questdlg( 'Merge these clusters?','Yes','No' );
        switch merge
            case 'Yes'
                bothLabels = [label1,label2];
                smallerLabel = min( bothLabels ); 
                if smallerLabel == bothLabels(1)
                    handles.labels(pts2) = label1;
                else
                    handles.labels(pts1) = label2;
                end
        end
        
        % update
        handles = st_clearSelectedData( handles );
        guidata( hObject,handles );

        
    % --- Executes on button press in deleteLabel.
    function deleteLabel_Callback(hObject, eventdata, handles)
        % deletes the clustering result, removing all previous labels

        disp( 'Click on a cluster to remove the label or data points' );
        axes( handles.mapplot );

        % get the nearest point to the selection
        pts = st_selectDataCluster( handles.mapplot,handles.labels );
        handles.selectedPoints(pts) = true;
        handles.selectedPointPlot = st_plotSelectedData( handles );
        deleteClust = questdlg('What do you want to delete?','Please choose an anser',...
            'Labels only','Labels & data points','Cancel','Cancel' );
        switch deleteClust
            case 'Labels only'
                handles.labels(pts) = 0;                
            case 'Labels & data points'
                handles = st_clearSelectedData( handles );
                handles.mapplot.Children.XData(pts) = [];
                handles.mapplot.Children.YData(pts) = [];
                handles.mapplot.Children.CData(pts,:) = [];
                handles = st_rmdata( handles,pts );
        end
        
        % update
        handles = st_clearSelectedData( handles );
        guidata( hObject,handles );
       
        
    % --- Executes on button press in finishManual.
    function finishManual_Callback(hObject, eventdata, handles)
        % hObject    handle to finishManual (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)

        % get all new clusters and the points belonging to them
        if sum( ~cellfun( @isempty,handles.manualClust ) > 0 )
            for i = 1:numel( handles.manualClust ) 
                if isvalid( handles.manualClust{i} )
                    position = handles.manualClust{i}.getPosition;
                    inpts = inpolygon( handles.projection(:,handles.plotDims(1)),...
                                handles.projection(:,handles.plotDims(2)),...
                                position(:,1),position(:,2) );

                    % create a new label for each new additional label
                    handles.labels(inpts) = max( handles.labels )+1;
                end
                delete( handles.manualClust{i} );
            end
        end
        
        % now hide the extra buttons    
        set( handles.newLabel,'enable','off' );
        set( handles.mergeLabel,'enable','off' );
        set( handles.deleteLabel,'enable','off' );
        set( handles.finishManual,'enable','off' );
        set( handles.autosort,'enable','on' );
        set( handles.projectdata,'enable','on' );
        set( handles.loaddata,'enable','on' );
        set( handles.replogt,'enable','on' );
        
        % update
        handles = update_labels( handles );
        handles.plotcolor = update_plot_colors( handles );
        st_update_plots( handles );
        guidata( hObject,handles );
    
        
    % --- Update all the labels
    function handles = update_labels( handles )
        % updates the labels so that new labels are added for each new clust
        
        % compute the means/covariances of the clusters
        proj = gather( handles.projection );
        [mu,sigma] = get_cluster_description( proj,handles.labels );
        
        % get the probabilities of the points belonging to the clusters,
        % assuming a gaussian mixture model
        [probabilities,priors] = compute_cluster_probabilities( proj,handles.labels,mu,sigma );
        
        % store the model and probabilities. Even if only a subset of
        % points are re-sorted, the probabilities for all points and
        % clusters will be recomputed. This provides a working model for
        % future spike sorting routines
        model = struct;
        model.mu = mu;
        model.Sigma = sigma;
        model.w = priors;
        handles.model.sortModel = model;
        handles.probabilities = probabilities;
                    
        
% ------------------------------------------------------------------------- 
% random callbacks/fcns not assigned
% -------------------------------------------------------------------------
    function projectdata_DeleteFcn(hObject,eventdata,handles)
        % hObject
        % eventdata
        % handles
        
        
    function Untitled_1_Callback(hObject, eventdata, handles)
        % hObject    handle to Untitled_1 (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)