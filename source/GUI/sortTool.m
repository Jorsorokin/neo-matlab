function varargout = sortTool( varargin )
    % code for sortTool.fig
    %
    % GUI for facilitating spike sorting of large databases. Inputs are optional,
    % as one can load in the various files through the GUI itself. This GUI
    % works equally well for any data that needs sorting (genomic, consumer preferences,
    % network graphs, etc.). The only difference is that some features may
    % be irrelevant (i.e. ISI and rasters, for instance).
    %
    % The GUI projects high-dimensional data using a variety of projection
    % methods that the user can specify. These projection methods are a
    % mixture of those written by the "FastICA" team, Laurens van der
    % Maaten (his awesome dimensionality redux. toolbox), and some written by
    % myself. The user can choose how many dimensions to project onto, 
    % as well as which two dimensions (i.e. which columns of the projected
    % data matrix) to visualize. Clicking "re-plot" will plot the
    % projections using a new viewpoint if specified.
    %
    % Following projections, one can choose to automatically sort
    % the clusters - and can further specify the number of clusters to sort
    % or allow an algorithm to solve for this number - using any of the following
    % clustering algorithms:
    % 
    % parametric:
    %   - Expectation Maximization for Gaussian Mixture Models (EM-GMM)
    %   - Expectation Maximization for T-distribution Mixture Models
    %     (EM-TMM...currently under development)
    %   - K-means (Km) 
    %   - Variational Bayes for GMMs (VB)
    % 
    % non-parametric:
    %   - Density-based spatial clustering for applications with noise (DBSCAN)
    %   - Hierarchical DBSCAN (HDBSCAN)
    %   - Spectral clustering
    %
    % Finally, one can also manually delete, merge, and create new clusters 
    % by clicking on the "Cut clusters" button, and can select individual 
    % data points or individual clusters for viewing associated waveforms, 
    % ISI, and raster plots by clicking the "Select points" and the
    % "Select cluster" buttons, respectively. To remove all cluster labels, 
    % click "Delete labels", and to close the GUI and extract the final
    % cluster labels, click "Finish". 
    %
    % INPUTS:
    %   ( data ) - the raw waveforms in column-order (each column is one
    %              observation, rows are variables). Can be 2- or
    %              3-dimensional, in which case the third dimension may be
    %              interpreted as channels
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
    %   ( R ) - a structure containing information on the projection / sorting: 
    %                   mapping         - the projection mapping structure (if applicable)
    %                   projMethod      - the projection method (PCA, ICA, LE,...)
    %                   sortMethod      - the sorting method (Km, EM-GMM, DBSCAN,... )
    %                   sortModel       - the sorting model structure (if applicable)
    %                   probabilities   - the probabilities of each spike
    %                                       belonging to each kth cluster
    %
    % Copywrite: Jordan Sorokin (Jorsorokin@gmail.com)
    % Last Edited: 6/22/2018
    
    %% =============== Figure Initiation ==================

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

    function sortTool_OpeningFcn(hObject, eventdata, handles, varargin)
        
        % Choose default command line output for sortTool
        handles.output = hObject;

        % set up the figure
        set( handles.figure1,'Name','Dimensionality Reduction & Sorting Tool',...
            'GraphicsSmoothing','off',...
            'CloseRequestFcn',@CloseMenuItem_Callback,...
            'NumberTitle','off',...
            'Visible','off',...
            'Toolbar','figure',...
            'MenuBar','none',...
            'WindowStyle','Normal' );
       
        uimenu( handles.figure1,...                     % creates a 'load' menu button
            'Label','Load',...
            'CreateFcn',@loaddata_OpeningFcn );
        
        uimenu( handles.figure1,...                     % creates "options" menu figure
            'Label','Options',...
            'Callback',@(hObject,eventdata)...
            sortTool('options_Callback',hObject,eventdata,guidata(hObject)) );
        
        % set up data parameters
        handles.data = nan;                             % the waveforms
        handles.availableData = false;                  % to avoid "gathering" gpuArray for each logical check
        handles.times = nan;                            % for plotting rasters/ISI
        handles.projection = nan;                       % the projection using the "projectMethod"    
        handles.labels = nan;                           % previous sorting labels
        handles.trials = nan;                           % n-dim vector used for plotting rasters
        handles.mask = nan;                             % mask vector (see "double_flood_fill.m")
        handles.location = nan;                         % vector specifying (x,y) location of each electrode
        handles.nLocationDim = 0;                       % the number of location dimensions
        handles.selectedPoints = nan;                   % points highlighted on the projection plot
        handles.selectedClusts = nan;                   % used for getting colors
        handles.manualClust = {};                       % will be added to while with manual cluster assignment
        handles.clusterIDs = uint8( 0 );
        
        % set up projection parameters
        handles.availableProjection = false;            % to avoid "gathering" gpuArray for each logical check
        handles.rotationDim = 1;                        % the dimension to rotate over
        handles.rotationAngle = 0.1;                    % radians to rotate by
        handles.arrows = [1 0; -1 0; 0 1; 0 -1];        % arrow positions for projection navigator
        handles.axis = [1 1 2 2];                       % axis (horizontal, vertical) for navigator
        handles.directions = [1 -1 1 -1];               % directions (left,right,up,down) for navigator
        handles.concatChans = false;                    % for concatenating channels or not for dim redux
        
        % set up sorting parameters
        handles.sortOptions.searchForK = 0;             % cluster # searching
        handles.sortOptions.rejectProb = 0.1;           % rejection probability
        handles.sortOptions.clustErr = 0;               % if 1, uses a metric to calculate clustering error based on densities
        handles.sortOptions.neighborType = 'kNN';       % neighborhood for Spectral/DBSCAN
        handles.sortOptions.laplaceNorm = 1;            % normalizes laplacian via random-walk
        handles.sortOptions.weightGraph = 0;            % RBF weighting of affinity matrix 
        handles.sortOptions.nNeighbors = 15;            % # neighbors for neighborhood graphs (including projections)
        handles.sortOptions.epsilon = 0.1;              % epsilon-radius neighborhood 
        handles.sortOptions.sigmaRBF = 1;               % SD of gaussian RBF kernel
        handles.sortOptions.K = 2;                      % defaults to 2 clusters
        handles.sortOptions.minpts = 5;                 % 5 pts min for density-based clustering
        handles.sortOptions.minclustsize = 5;           % any clusters with < this are set to 0
        handles.sortOptions.outlierThresh = 0.9;        % for HDBSCAN
        handles.sortOptions.clusterMetric = 'Silhouette'; % For measuring cluster quality
        handles.sortOptions.tryMerge = false;           % tries to merge clusters upon "refine clusts" pushbutton if true

        % create R structure for ouput
        handles.R = struct; 
        handles.R.projectMethod = 'PCA';                % default projection method
        handles.R.sortMethod = 'EM-GMM';                % default sorting method
        handles.R.mapping = nan;                        % the projection mapping
        handles.R.keptPts = nan;                        % updates if we delete points
        handles.R.clusterQuality = nan;                 % stores the quality of the clustering 
        
        % plotting colors for different labels
        handles.allPlotColor = [ [.6, .6, .6];...       % gray              (0)
                                 [0.1, 0.74, 0.95];...  % deep sky-blue     (1)
                                 [0.95, 0.88, 0.05];... % gold/yellow       (2)
                                 [0.80, 0.05, 0.78];... % magenta           (3)
                                 [0.35, 0.8, 0.3];...   % lime green        (4)
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
                                 [0.53, 0.85, 0.98];... % light sky blue    (15)
                                 [0, 0.6, 0.1];...      % forest green      (16)
                                 [0.65, 0.9, 0.6];...   % light green       (17)
                                 [0.85, 0.6, 0.88];...  % light purple      (18)
                                 [0.90, 0.7, 0.7];...   % light red         (19)
                                 [0.25, 0.25, 0.5];...  % dark blue         (20)
                                ];
        handles.plotcolor = nan;    
        
        % GET THE OPTIONAL INPUTS
        p = st_check_inputs( varargin );
        
        handles.data = p.data;
        if any( ~isnan( gather( handles.data ) ) )
            handles.availableData = true;
            nPoints = size( handles.data,2 );
        end
        
        handles.projection = p.projection;
        if any( ~isnan( gather( handles.projection ) ) )
            handles.availableProjection = true;
            
            % create the plotting matrices
            [nPoints,handles.nDim] = size( handles.projection );
            handles.plotDims = zeros( handles.nDim,2 );
            handles.plotDims(1,1) = 1; handles.plotDims(2,2) = 1;
            handles.loadingMatrix = handles.plotDims;
        end
        
        if any( [handles.availableData,handles.availableProjection] )
            handles.R.keptPts = true( 1,nPoints );
        end
        
        handles.labels = p.labels;
        handles.times = p.times;
        handles.trials = p.trials;
        handles.mask = p.mask;
        handles.location = p.location;
        if ~all( isnan( handles.location ) )
            handles.nLocationDim = size( handles.location,2 );
        end
        
        if all( isnan( handles.trials ) ) && ~any( isnan( handles.times ) )
            handles.trials = ones( 1,numel( handles.times ),class( handles.times ) );
        end
        
        clear p
               
        % get the selectedPoints if data available
        if ~isnan( handles.labels )
            handles.selectedPoints = false( 1,numel( handles.labels ) );
            handles.labels = uint8( handles.labels );
            handles.clusterIDs = unique( handles.labels );
            handles.selectedClusts = false( 1,numel( handles.clusterIDs ) );
            if any( isnan( handles.R.keptPts ) )
                handles.R.keptPts = true( 1,numel( handles.labels ) );
            end
        end
        
        % INITIATE THE PLOTS
        set( handles.mapplot,'tickdir','out','box','off','fontsize',8,...
            'xcolor','k','ycolor','k','NextPlot','add');
        handles.waveformplot = nan;
        handles.isiplot = nan;
        handles.rasterplot = nan;
        handles.xcorrplot = nan;
        handles.driftplot = nan;
        handles.qualityplot = nan;
        handles.dvdtplot = nan;
        handles.loadingplot = [];
        handles.plotcolor = update_plot_colors( handles );
        
        % plot the projections if provided
        if handles.availableProjection
            st_plot_projections( handles );
        end      
        
        % turn on the main figure
        handles.figure1.Visible = 'on';
        
        % Update handles structure 
        %setappdata( 0,'HandleMainGUI',hObject );
        guidata( hObject, handles );
        
        % keep track of the GUI visibility before running the OutputFcn
        waitfor( handles.figure1,'Visible','off' );
        

    function varargout = sortTool_OutputFcn( hObject, eventdata, handles )
        varargout{1} = handles.labels;
        varargout{2} = handles.projection;
        varargout{3} = handles.R;
        
        close all
%         if isfield( handles,'optionsFig' )
%             close( handles.optionsFig );
%         end
%         
%         if ~isa( handles.waveformplot,'double' ) && isvalid( handles.waveformplot )
%             close( handles.waveformplot );
%         end
%         if ~isa( handles.isiplot,'double' ) && isvalid( handles.isiplot )
%             close( handles.isiplot );
%         end
%         if ~isa( handles.rasterplot,'double' ) && isvalid( handles.rasterplot )
%             close( handles.rasterplot.Parent );
%         end
%         if ~isa( handles.qualityplot,'double' ) && isvalid( handles.qualityplot )
%             close( handles.qualityplot.Parent );
%         end
%         
%         close( handles.figure1 );
                 

    function CloseMenuItem_Callback(hObject, eventdata)
        handles = guidata( hObject );
        delete( handles.figure1 );
        if isfield( handles,'optionsFig' )
            delete( handles.optionsFig );
        end
        close all
        

    function savebutton_Callback(hObject, eventdata, handles)
        % closes the GUI   
        handles.figure1.Visible = 'off';

    
    function options_Callback( hObject,eventdata,handles )
        
        % create new optionsFig if not available
        if ~isfield( handles,'optionsFig' )
            handles.optionsFig = st_create_optionsFig;
            handles.optionsFig.CloseRequestFcn = @options_CloseFcn;
        end
        
        % open the optionsFig
        handles.optionsFig.Visible = 'on';

        guidata( hObject,handles );


    function options_CloseFcn( hObject,eventdata )
        
        % retrieve handles
        mainFig = hObject.Parent.Children(2);
        handles = guidata( mainFig ); % retrieves correct handles object

        % make options invisible
        handles.optionsFig.Visible = 'off';

        % store the values into the relevant variables in 'handles'
        h = findobj( handles.optionsFig,'Tag','General' );
        handles.sortOptions.searchForK = get( findobj( h,'Tag','searchForK' ),'Value' );
        handles.sortOptions.rejectProb = str2double( get( findobj( h,'Tag','probBox' ),'String' ) );
        handles.sortOptions.K = str2double( get( findobj( h,'Tag','kBox' ),'String' ) );
        handles.sortOptions.outlierThresh = str2double( get( findobj( h,'Tag','outlierBox' ),'String' ) );
        handles.sortOptions.tryMerge = get( findobj( h,'Tag','tryMerge' ),'Value' );
        clusterMetric = findobj( h,'Tag','clusterMetricType' );
        clusterMetric = clusterMetric.SelectedObject.String;
        if ~strcmp( handles.sortOptions.clusterMetric,clusterMetric )
            handles.sortOptions.clusterMetric = clusterMetric;
            if ~isa( handles.qualityplot,'double' ) && isvalid( handles.qualityplot )
                cla( handles.qualityplot.Children(1) );
                cla( handles.qualityplot.Children(2) );
            end
        end

        h = findobj( handles.optionsFig,'Tag','Graphs' );
        neighborType = findobj( h,'Tag','neighborType' );
        handles.sortOptions.neighborType = neighborType.SelectedObject.String;
        handles.sortOptions.laplaceNorm = get( findobj( h,'Tag','laplaceNorm' ),'Value' );
        handles.sortOptions.weightGraph = get( findobj( h,'Tag','weightGraph' ),'Value' );
        handles.sortOptions.nNeighbors = str2double( get( findobj( h,'Tag','nnBox' ),'String' ) );
        handles.sortOptions.epsilon = str2double( get( findobj( h,'Tag','epsBox' ),'String' ) );    
        handles.sortOptions.sigmaRBF = str2double( get( findobj( h,'Tag','sigmaBox' ),'String' ) ); 
        handles.sortOptions.minclustsize = str2double( get( findobj( h,'Tag','minclustBox' ),'String' ) );
        
        guidata( mainFig,handles );
        

    function keyPress_Callback( hObject,eventdata,handles )        
        switch eventdata.Key
            case 'leftarrow'
                selectPrevClust( hObject,handles );
            case 'rightarrow'
                selectNextClust( hObject,handles );
        end

    %% ================ loading & projections ================
    
    function loaddata_OpeningFcn( hObject,eventdata )
        handles = guidata( hObject );
        loadNames = {'data','times','projection','labels','trials','mask','location'};
        for j = 1:numel( loadNames )
            uimenu( hObject,... % creats submenu variables to load
                'Label',loadNames{j},...
                'UserData',handles,...
                'Callback',@loaddata_Callback );
        end
        
        % update the handles object
        guidata( hObject,handles );
    
        
    function loaddata_Callback(hObject, eventdata)
        
        % retrieve most recent handles object, clear plots
        handles = guidata( hObject );
        st_clear_plots( handles );

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
            x = evalin( 'base',name ); % loads in data into 'x'
            switch hObject.Label
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
            'LDA (linear discriminant analysis',...
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


    function mappingmethod_Callback(hObject, eventdata, handles)
        % selects a new mapping method

        allmethods = {'PCA','KernelPCA','ProbPCA','NCA','LPP','LDA','ICA','LLE',...
                        'HessianLLE','NPE','SPE','MVU','FastMVU','Laplacian',...
                        'DiffusionMaps','SNE','tSNE','Isomap','Sammon','Autoencoder'};
        value = get( hObject,'value' );
        handles.R.projectMethod = allmethods{value};

        % update the GUI
        guidata( hObject, handles );

        
    function projectData_Callback(hObject,eventdata,handles)
        % projects waveforms using the chosen mapping method

        % do the projections
        nDim = findobj( hObject.Parent,'tag','nDim' );
        nDim = str2double( nDim.String );
        if isnan( nDim )
            fprintf( 'Please specify # of dimensions to project onto\n' );
            return
        elseif nDim > size( handles.data,2 )
            fprintf( '# of Dimensions must be < # of data points\n' );
            return
        else
            handles.nDim = nDim;
        end

        [handles.projection,handles.R.mapping] = st_project_data( handles );
        handles.nDim = size( handles.projection,2 ); 
        handles.availableProjection = true;
        
        % update the projection viewer based on the # of
        % dimensions that have been kept
        totalDims = handles.nDim + handles.nLocationDim;
        if totalDims == 1
            handles.plotDims = 1;
        else
            handles.plotDims = zeros( totalDims,2 );
            handles.plotDims(1,1) = 1;
            handles.plotDims(2,2) = 1;
        end
        
        handles.loadingMatrix = handles.plotDims;
        textObj = findobj( handles.figure1.Children,'Tag','rotationDimText' );
        textObj.String = 'dim 1';
        handles.rotationDim = 1;
        
        % now plot onto the main figure
        st_plot_projections( handles );
        st_plotSelectedData( handles );
        st_update_loadingplot( handles,1 );
        st_update_loadingplot( handles,2 );

        % update the GUI
        guidata( hObject, handles );
        

    function nDim_Callback(hObject, eventdata, handles)
        % specifies the number of dimensions to project onto using the specified mapping method

        handles.nDim = str2double( get( hObject,'String' ) );
        guidata( hObject,handles );    
        
    
    function concatChans_Callback( hObject,eventdata,handles )
        % sets the flag to concatenates channels before dim redux or not
        
        handles.concatChans = hObject.Value;
        guidata( hObject,handles )
    
        
    %% ================ plotting procedures ================
    
    function projectionNavigator_createFcn( hObject,eventdata,handles )
        % creates the projection navigator window / graphics
        
        % plot a compass onto the axis
        hold on
        plotHandle = quiver( hObject,...
            [0;0;0;0],[0;0;0;0],[0.5;-0.5;0;0],[0;0;0.5;-0.5],...
            'color','w' );
        plotHandle.PickableParts = 'none';
        plotHandle.HitTest = 'off';     
        
        hObject.XLim = [-0.5,0.5];
        hObject.YLim = [-0.5,0.5];
        
        %handles.projectionNavigator = hObject;
        guidata( hObject,handles );
        
            
    function projectionNavigator_clickFcn( hObject,eventdata,handles )
        % allows the user to rotate through different views of the
        % projections by changing the loadings for different projection
        % dimensions in the loading matrix
        
        if ~handles.availableProjection
            return
        end
        
        % check where the user clicked
        mousePos = eventdata.IntersectionPoint(1:2);
        [~,closestArrow] = min( compute_pairwise_dist( mousePos,handles.arrows ) );
        
        % change the value of the handles.plotDim variable depending 
        % on which dimension we are rotating over, and which direction
        % (left/right vs. up/down)   
        loadingMatrix = handles.loadingMatrix;
        rotationAxis = handles.axis(closestArrow);
        rotationValue = loadingMatrix(handles.rotationDim,rotationAxis);
        rotationAmnt = handles.rotationAngle * handles.directions(closestArrow);
        rotationValue = rotationValue + rotationAmnt;
        handles.loadingMatrix(handles.rotationDim,rotationAxis) = rotationValue;
        handles.plotDims = sin( handles.loadingMatrix );
        
        % replot the scatter and update the loading plot 
        st_update_loadingplot( handles,rotationAxis );
        
        if rotationAxis == 1
            handles.mapplot.Children.XData = handles.projection * handles.plotDims(:,1);
        else
            handles.mapplot.Children.YData = handles.projection * handles.plotDims(:,2);
        end
        
        % update handles
        guidata( hObject,handles );
     
        
    function prevDim_callback( hObject,eventdata,handles )
        % moves the projection navigator to the next dimension of the
        % projections for rotating the scatter
        
        handles.rotationDim = max( handles.rotationDim - 1,1 );        
        textObj = findobj( handles.figure1.Children,'Tag','rotationDimText' );
        textObj.String = ['dim ',num2str( handles.rotationDim )];
        guidata( hObject,handles );

            
    function nextDim_callback( hObject,eventdata,handles )
        % moves the projection navigator to the previous dimension of the
        % projections for rotating the scatter
        
        handles.rotationDim = min( handles.rotationDim + 1,handles.nDim + handles.nLocationDim );
        textObj = findobj( handles.figure1.Children,'Tag','rotationDimText' );
        textObj.String = ['dim ',num2str( handles.rotationDim )];
        guidata( hObject,handles );
    
       
    function viewProjectionLoading_callback( hObject,eventdata,handles )
        % plots the loading weights of the dimensions of the projections 
        % for each of the two axes of the scatter plot
        
        if ~handles.availableProjection
            return
        end
        
        if ~isempty( handles.loadingplot ) && handles.loadingplot.isvalid 
            return
        end
        
        % create the figure
        handles.loadingplot = st_create_loadingplot( handles );
        guidata( hObject,handles );
        
        
    function resetProjectionLoading_callback( hObject,eventdata,handles )
        % resets the loadings so that the first and second dimensions are 
        % plotted completely onto the x-axis / y-axis 
        
        handles.plotDims = 0*handles.plotDims;
        handles.plotDims(1,1) = 1;
        handles.plotDims(2,2) = 1;
        handles.loadingMatrix = handles.plotDims;
        
        % update the scatter
        st_update_scatter( handles );
        st_update_loadingplot( handles,1 );
        st_update_loadingplot( handles,2 );
        guidata( hObject,handles );
        
        
    function randomizeProjections_callback( hObject,eventdata,handles )
        % projects the high dimensional data onto the two dimensional
        % figure using a random projection of the features
        
        nDims = size( handles.plotDims,1 );
        handles.plotDims = rand( nDims,2 )*2 - 1;
        handles.loadingMatrix = handles.plotDims;
        
        % update the scatter
        st_update_scatter( handles )
        st_update_loadingplot( handles,1 );
        st_update_loadingplot( handles,2 );
        guidata( hObject,handles );
                
        
    function selectdata_Callback(hObject, eventdata, handles)
        % allows user to select individual data points in the main plot  

        if ~handles.availableProjection
            return
        end

        % get the selection value, plot the selected points
        pts = st_selectDataLasso( handles.mapplot );
        if iscell( pts )
            pts = pts{end};
        end
        
        if isempty( pts )
            return
        end
        
        % plot the selected points
        handles.selectedPoints(pts) = true;
        handles.selectedClusts(ismember( handles.clusterIDs,handles.labels(pts) )) = true;
        st_plotSelectedData( handles );
        
        % update 
        st_update_plots( handles );
        guidata( hObject,handles );
        
    
    function selectclust_Callback(hObject,eventdata,handles)
        % allows user to select individual clusters in the main plot

        if ~handles.availableProjection
            return
        end
        
        % get the selected point and those belonging to same cluster
        pts = st_selectDataCluster( handles.mapplot,handles.labels );
        if iscell( pts )
            pts = pts{end};
        end
        
        if isempty( pts )
            return
        end
        
        % plot selected cluster
        handles.selectedPoints(pts) = true;
        handles.selectedClusts(ismember( handles.clusterIDs,handles.labels(pts(1)) )) = true;
        st_plotSelectedData( handles );
        
        % update
        st_update_plots( handles );
        guidata( hObject,handles );
        
        
    function selectNextClust( hObject,handles )
        % allows user to slect previous (sequentially) cluster
               
        % check if points selected are all in one cluster and entire
        % cluster is selected
        allID = handles.clusterIDs;
        uID = allID( handles.selectedClusts );
        if numel( uID ) ~= 1 || any( handles.labels(~handles.selectedPoints) == uID )
            return
        end    

        % deselect all points
        st_clearSelectedData( handles );
        handles.selectedPoints(handles.selectedPoints) = false;
        
        % select the next cluster
        clustInd = strfind( allID,uID );
        if clustInd == numel( allID )
            clustInd = 0;
        end
        nextClust = allID(clustInd + 1);
        handles.selectedPoints(handles.labels==nextClust) = true;
        handles.selectedClusts(handles.selectedClusts) = false;
        handles.selectedClusts(clustInd+1) = true;
        
        % update
        st_plotSelectedData( handles );
        st_update_plots( handles );
        guidata( hObject,handles );
        
        
    function selectPrevClust( hObject,handles )
        % allows user to slect previous (sequentially) cluster
        
        % check if points selected are all in one cluster and entire
        % cluster is selected
        allID = handles.clusterIDs;
        uID = allID( handles.selectedClusts );
        if numel( uID ) ~= 1 || any( handles.labels(~handles.selectedPoints) == uID )
            return
        end    
        
        % deselect all points
        st_clearSelectedData( handles );
        handles.selectedPoints(handles.selectedPoints) = false;
        
        % select the previous cluster
        clustInd = strfind( allID,uID );
        if clustInd == 1
            clustInd = numel( allID ) + 1;
        end
        nextClust = allID(clustInd - 1);
        handles.selectedPoints(handles.labels==nextClust) = true;       
        handles.selectedClusts(handles.selectedClusts) = false;
        handles.selectedClusts(clustInd - 1) = true;
        
        % update
        st_plotSelectedData( handles );
        st_update_plots( handles );
        guidata( hObject,handles );
        
        
    function cleardata_Callback(hObject, eventdata, handles)

        % clear all data selection
        st_clearSelectedData( handles );
        handles.selectedPoints(handles.selectedPoints) = false;
        handles.selectedClusts(handles.selectedClusts) = false;
        
        % update
        guidata( hObject,handles );
        
                    
    function colors = update_plot_colors( handles )
        
        % check if any labels (i.e. any data)
        if any( isnan( handles.labels ) )
            colors = nan;
            return
        else
            colors = zeros( numel( handles.labels ),3,'single' );
        end
        
        % update plotting colors for each data point
        nColors = size( handles.allPlotColor,1 ) - 1;
        for i = handles.clusterIDs
            idx = (handles.labels == i);
            thisColor = min( i,mod( i,nColors ) + uint8(i ~= 0) );
            colors(idx,:) = repmat( single( handles.allPlotColor(thisColor+1,:) ),sum( idx ),1 ); % since uID==0 is first color
        end
            
               
    function plotspikes_Callback(hObject, eventdata, handles)
        
        % create a new figure if needed
        if isa( handles.waveformplot,'double' ) || ~isvalid( handles.waveformplot )
            handles.waveformplot = st_create_waveplot( handles ); 
        end
            
        % plot
        switch handles.waveformplot.Visible
            case 'on'
                arrayfun( @cla,handles.waveformplot.Children );
                handles.waveformplot.Visible = 'off'; % hide the plot
            otherwise
                st_plot_waveforms( handles )
                handles.waveformplot.Visible = 'on';
        end
        
        % update
        guidata( hObject,handles );    
        
        
    function plotisi_Callback(hObject, eventdata, handles)
        
        % create a new figure if needed
        if isa( handles.isiplot,'double') || ~isvalid( handles.isiplot )
            handles.isiplot = st_create_isiplot();
        end
        
        % check if visible; if not, plot
        switch handles.isiplot.Visible
            case 'on'
                cla( handles.isiplot.Children );
                handles.isiplot.Visible = 'off';
            otherwise
                st_plot_isi( handles );
                handles.isiplot.Visible = 'on';
        end
        
        % update
        guidata( hObject,handles );  
        
        
    function plotraster_Callback(hObject, eventdata, handles)
        
        % create a new figure if needed
        if isa( handles.rasterplot,'double' ) || ~isvalid( handles.rasterplot )
            handles.rasterplot = st_create_rasterplot( handles );
        end
        
        % check if Visible; if not, plot the raster
        switch handles.rasterplot.Visible
            case 'on'
                arrayfun( @cla,handles.rasterplot.Children );
                handles.rasterplot.Visible = 'off';
            otherwise
                st_plot_raster( handles );
                handles.rasterplot.Visible = 'on';
        end
        
        % update
        guidata( hObject,handles );  

        
    function plotxcorr_Callback(hObject, eventdata, handles)
        
        % create a new figure if needed
        if isa( handles.xcorrplot,'double' ) || ~isvalid( handles.xcorrplot )
            handles.xcorrplot = st_create_correlogramplot();
        end
        
        % check if Visible; if not, plot the cross-correlograms
        switch handles.xcorrplot.Visible
            case 'on'
                handles.xcorrplot.Visible = 'off';
            otherwise
                delete( handles.xcorrplot.Children );
                st_plot_correlograms( handles );
                handles.xcorrplot.Visible = 'on';
        end
        
        % update
        guidata( hObject,handles );  
        
        
    function plotdrift_Callback(hObject, eventdata, handles)
        
        % create a new figure if needed
        if isa( handles.driftplot,'double' ) || ~isvalid( handles.driftplot )
            handles.driftplot = st_create_driftplot( handles );
        end
        
        % check if Visible; if not, plot the cross-correlograms
        switch handles.driftplot.Visible
            case 'on'
                handles.driftplot.Visible = 'off';
            otherwise
                arrayfun( @cla,handles.driftplot.Children );
                st_plot_drift( handles );
                handles.driftplot.Visible = 'on';
        end
        
        % update
        guidata( hObject,handles );
        

    function plotdvdt_Callback(hObject, eventdata, handles)
        
        % create a new figure if needed
        if isa( handles.dvdtplot,'double' ) || ~isvalid( handles.dvdtplot )
            handles.dvdtplot = st_create_dvdtplot( handles );
        end
        
        % check if Visible; if not, plot the cross-correlograms
        switch handles.dvdtplot.Visible
            case 'on'
                handles.dvdtplot.Visible = 'off';
            otherwise
                cla( handles.dvdtplot.Children );
                st_plot_dvdt( handles );
                handles.dvdtplot.Visible = 'on';
        end
        
        % update
        guidata( hObject,handles );    
        
        
    function plotClusterQuality_Callback( hObject,eventdata,handles )

        % check if cluster quality is available
        if isa( handles.qualityplot,'double' ) || ~handles.qualityplot.isvalid
            handles.qualityplot = st_create_qualityPlot( handles );        
        end
        
        % plot the cluster quality for the clusters
        switch handles.qualityplot.Visible 
            case 'on'
                handles.qualityplot.Visible = 'off';
            case 'off'
                st_plot_clusterQuality( handles );
        end
        
        % update
        guidata( hObject,handles );
        
                
    function mapplot_mouseclick( hObject,eventdata,handles )
        % allows the user to right-click on the scatter plot and delete
        % points/labels for each point separately
        
        % check if any data projected
        if isempty( hObject.Children )
            return
        end
        
        % check for right-click only
        if ~strcmp( hObject.Parent.SelectionType,'alt' )
            return
        end
        XData = hObject.Children.XData;
        YData = hObject.Children.YData;
        xy = get( hObject,'CurrentPoint' );
        xy = xy(1,1:2);
        ax = axis;
        dx = ax(2) - ax(1);
        dy = ax(4) - ax(3);
        [~,closestPt] = min( sqrt( ((XData - xy(1))/dx).^2 + ((YData - xy(2))/dy).^2 ) );

        % draw a circle around the point
        pos = [XData(closestPt),YData(closestPt)];
        hold on
        circle = plot( hObject,pos(1),pos(2),'wo','linewidth',2 );

        % see if user wants to delete point or label
        answer = questdlg( 'Delete point?',' ',...
            'Label only','Label & point','Cancel','Cancel' );

        % remove the circle
        delete( circle );      
        hold off        

        % remove label and/or point      
        switch answer
            case 'Label only'
                handles.labels(closestPt) = 0;

            case 'Label & point'
                warning 'off' 
                
                handles.data(:,closestPt,:) = [];
                handles.labels(closestPt) = [];
                hObject.Children.SizeData(closestPt) = [];
                hObject.Children.XData(closestPt) = [];
                hObject.Children.YData(closestPt) = [];
                hObject.Children.CData(closestPt,:) = [];
                handles = st_rmdata( handles,closestPt );
                
                warning 'on'
        end

        % update the handles
        handles = update_sortmodel( handles );
        handles.plotcolor = update_plot_colors( handles );
        st_update_plots( handles );
        guidata( hObject,handles ); 


    function deletebutton_Callback(hObject,eventdata,handles)
        if isnan( handles.labels )
            return
        end
        
        % clear previous labels of highlighted pts, update the plots
        if ~any( handles.selectedPoints )
            pts = 1:size( handles.projection,1 );
        else
            pts = handles.selectedPoints;
        end
        handles.labels(pts) = 0;
        handles.clusterIDs = unique( handles.labels );
        handles.selectedClusts = false( 1,numel( handles.clusterIDs ) );
        handles.plotcolor = update_plot_colors( handles );
        st_update_plots( handles );
        guidata( hObject,handles );        

    %% ================== Sorting =====================

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
            'HDBSCAN',...
            'Spectral clustering'} );

        % update the GUI
        guidata( hObject, handles ); 
        
        
    function sortmethod_Callback(hObject,eventdata,handles)
        % select a new sorting method

        allmethods = {'EM-GMM','mEM-GMM','EM-TMM','Km',...
            'VB','mVB','DBSCAN','HDBSCAN','Spectral'};
        handles.R.sortMethod = allmethods{hObject.Value};

        % update the GUI
        guidata( hObject, handles );

        
    function autosort_Callback(hObject, eventdata, handles)        
        % automatically sort the data projections via the specified sort method

        warning off;
        
        % first project the data 
        if ~handles.availableData && ~handles.availableProjection
            disp( 'No data or projections available. Please load data!' );
            return
        end
        
        if ~handles.availableProjection
            disp( 'Please project data first' );
            return
        end

        % get the selected points, if any
        if any( handles.selectedPoints )
            pts = handles.selectedPoints;
            saveModel = true;
        else
            continueSorting = questdlg( 'Sort all points?','','Yes','No','Yes' );
            if ~strcmp( continueSorting,'Yes' )
                return
            end
            pts = ~handles.selectedPoints;
            saveModel = false;
        end
        
        % request for the appropriate cluster number for parametric
        % clustering algorithms
        switch handles.sortOptions.searchForK
            case 0
                K = handles.sortOptions.K; 
            case 1
                switch handles.R.sortMethod
                    case {'mEM-GMM','EM-GMM','EM-TMM','VB','mVB','Km'}
                        if ~handles.availableProjection
                            disp( 'Please project your data first or load projections!' );
                            return
                        end
                        [~,K] = findClustNum( handles.projection(pts,:)',2,50 );
                    otherwise
                        K = 1;                        
                end
        end
        
        % compute mask if necessary
        switch handles.R.sortMethod
            case {'mVB','mEM-GMM'}
                
                if any( isnan( handles.mask ) )
                    disp( 'must include a masking matrix for masked-EM and masked-VB sorting' );
                    return
                end

                % compute the projections
                projections = compute_spike_features( gather( permute( handles.data(:,pts,:),[2,1,3] ) ),...
                    handles.nDim,handles.R.projectMethod,handles.mask,[],false );

                % compute the expanded mask
                mask = zeros( size( projections ) );
                for c = 1:size( handles.mask,1 )
                    inds = c*handles.nDim-handles.nDim+1:c*handles.nDim;
                    mask(:,inds) = repmat( handles.mask(c,pts)',1,handles.nDim );
                end 

            otherwise                    
                projections = handles.projection(pts,:);
                mask = nan;
        end
        
        % now do the actual clustering using the method specified
        [labels,model] = sort_clusters( projections,K,handles.R.sortMethod,...
                                                        'mask',mask,...
                                                        'neighbors',handles.sortOptions.nNeighbors,...
                                                        'sigma',handles.sortOptions.sigmaRBF,...
                                                        'neighborType',handles.sortOptions.neighborType,...
                                                        'eps',handles.sortOptions.epsilon,...
                                                        'kernelWeighting',handles.sortOptions.weightGraph,...
                                                        'minclustsize',handles.sortOptions.minclustsize,...
                                                        'outlierThresh',handles.sortOptions.outlierThresh );
                                                    
        % change the labels to contain the original label plus consecutive
        % labels for each new cluster
        uID = handles.clusterIDs;
        if max( uID ) > 0
            if ~any( handles.selectedPoints )
                labels = labels + min( uID(uID > 0) ) - 1; % subtract 1 to start with the minimum of the previous label
            else 
                prevClusts = uID(ismember( uID,handles.labels(pts) ));
                newClusts = unique( labels );
                nnzClusts = newClusts(newClusts > 0);
                if any( newClusts > 0 )
                    nPrev = numel( prevClusts );
                    nNew = nnz( newClusts > 0 );
                    val = uint8( min( [nPrev,nNew] ) );
                    for j = 1:val
                        tempIDX = ismember( labels,nnzClusts(j) );
                        labels(tempIDX) = prevClusts(j);
                    end
                    
                    if nNew > val
                        tempIDX = ismember( labels,nnzClusts(val+1:end) );
                        labels(tempIDX) = labels(tempIDX) - 1 + max( uID );
                    end
                end
            end
        end
        
        % replace the model ONLY IF no individual points are selected
        if saveModel % all will be false if none are specifically selected
            handles.R.sortModel = model;
        end
        
        % update all probabilities and labels in the handles struct
        handles.labels(pts) = labels;
        handles.clusterIDs = unique( handles.labels );
        if any( handles.selectedClusts )
            handles.selectedClusts(ismember( handles.clusterIDs,unique( labels ) )) = true;
        end
        
        handles = update_sortmodel( handles );
        handles.plotcolor = update_plot_colors( handles );
        st_update_plots( handles );
        guidata( hObject,handles );
        
        warning on;
    
        
    function refineClusts_Callback( hObject,eventdata,handles )
        % dynamically split/merge clusters after sorting based on 
        % chi-squared tests for chi-squared distributions of each cluster

        % check if we've sorted
        if any( isnan( handles.labels ) ) || all( handles.labels == 0 )
            fprintf( 'Please sort clusters first\n' );
            return
        end

        % try splitting clusters first       
        fprintf( 'attempting to split clusters' );
        handles.labels = refine_clusters( handles.data,handles.labels );
%         uID = handles.clusterIDs(handles.clusterIDs > 0);
%         maxID = max( uID );
%         alpha = 0.01 / numel( uID ); % the corrected threshold for significance
%         for i = uID
%             fprintf( '.' );
%             inds = (handles.labels == i);
%             if nnz( inds ) <= handles.nDim
%                 continue
%             end
%             
%             [~,~,v] = svd( handles.projection(inds,:)','econ' );
%             v = v(:,1:min(size( handles.projection,2 ),3));
%             [newLabels,~,pval] = try_cluster_split( v );
%             if max( newLabels ) == 2 && pval < alpha
%                 maxID = maxID + 1;
%                 newInds = (newLabels == 2);
%                 newLabels(newInds) = maxID;
%                 newLabels(~newInds) = i;
%                 handles.labels(inds) = newLabels;
%             end
%         end
%         fprintf( '\n' );
%         
%         % try merging clusters if try merge is true
%         if handles.sortOptions.tryMerge
%             fprintf( 'attempting to merge clusters...\n' );
%             merge = try_cluster_merge( handles.projection,handles.labels );
%         end
        
        % update plots
        handles = update_sortmodel( handles );
        handles.clusterIDs = unique( handles.labels );
        handles.plotcolor = update_plot_colors( handles );
        st_update_plots( handles );
        guidata( hObject,handles );
        
        
    function manualsort_Callback(hObject, eventdata, handles)
        % initiates the functionality for manual sorting (cluster cutting)

        % clear any data selection
        st_clearSelectedData( handles );
        handles.selectedPoints(handles.selectedPoints) = false;
        
        % check if data has been plotted
        axes( handles.mapplot );
        if isempty( handles.mapplot.Children )
            if handles.availableProjection
                st_plot_projections( handles );
            end
        end        
        
        % set visibility of the new buttons to "on"
        set( handles.newLabel,'enable','on' );
        set( handles.mergeLabel,'enable','on' );
        set( handles.deleteLabel,'enable','on' );
        set( handles.finishManual,'enable','on' );
        set( handles.autosort,'enable','off' );
        set( handles.projectdata,'enable','off' );
        set( handles.refineClusts,'enable','off' );

        % update
        guidata( hObject,handles );
        
        
    function newLabel_Callback(hObject, eventdata, handles)
        % manually select a new label by lasso-ing around data points

        fprintf( 'Select a new cluster\n' );
        axes( handles.mapplot );

        % draw lasso
        clust = imfreehand( 'Closed',1 );
        nClust = sum( ~cellfun( @isempty,handles.manualClust ) );
        handles.manualClust{nClust+1} = clust;
        
        % update
        guidata( hObject,handles );

        
    function mergeLabel_Callback(hObject, eventdata, handles)
        % merges two clusters together 

        fprintf( 'Click on two clusters to merge\n' );
        axes( handles.mapplot );
        
        % get the points and associated labels with each cluster
        pts = [];
        label = inf;
        msg = msgbox( 'Click "OK" when you are finished selecting clusters to merge' );
        msg.Position(1) = msg.Parent.ScreenSize(3)*0.5;
        while true
            [newPts,newLabels] = st_selectDataCluster( handles.mapplot,handles.labels );
            
            % don't add these points if the msg box isn't valid
            if ~isvalid( msg )
                break
            end
            
            % add pts to running total
            pts = [pts,newPts];
            label = min( label,newLabels );
            handles.selectedPoints(pts) = true;
            st_plotSelectedData( handles );
        end
               
        % merge if selected "yes"
        merge = questdlg( 'Merge these clusters?','Yes','No' );
        switch merge
            case 'Yes'
                handles.labels(pts) = label;
        end
        
        % update
        st_clearSelectedData( handles );
        handles.selectedPoints(handles.selectedPoints) = false;
        guidata( hObject,handles );

        
    function deleteLabel_Callback(hObject, eventdata, handles)
        % deletes the clustering R, removing all previous labels

        disp( 'Click on a cluster to remove the label or data points' );
        axes( handles.mapplot );

        % get the nearest point to the selection
        pts = st_selectDataCluster( handles.mapplot,handles.labels );
        handles.selectedPoints(pts) = true;
        st_plotSelectedData( handles );
        deleteClust = questdlg('What do you want to delete?','Please choose an answer',...
            'Labels only','Labels & data points','Cancel','Cancel' );
        switch deleteClust
            case 'Labels only'
                handles.labels(pts) = 0;     
            case 'Labels & data points'
                warning 'off'
                
                st_clearSelectedData( handles );
                handles.selectedPoints(handles.selectedPoints) = false;
                handles.mapplot.Children.XData(pts) = [];
                handles.mapplot.Children.YData(pts) = [];
                handles.mapplot.Children.CData(pts,:) = [];
                handles.mapplot.Children.SizeData(pts) = [];
                handles = st_rmdata( handles,pts );
                
                warning 'on'
        end
        
        % update
        st_clearSelectedData( handles );
        handles.selectedPoints(handles.selectedPoints) = false;
        guidata( hObject,handles );
       
        
    function finishManual_Callback(hObject, eventdata, handles)

        % get all new clusters and the points belonging to them
        if sum( ~cellfun( @isempty,handles.manualClust ) > 0 )
            scatterPoints = [handles.mapplot.Children(end).XData; handles.mapplot.Children(end).YData];
            maxID = handles.clusterIDs(end);
            for i = 1:numel( handles.manualClust ) 
                if isvalid( handles.manualClust{i} )
                    position = handles.manualClust{i}.getPosition;
                    inpts = inpolygon( scatterPoints(1,:),scatterPoints(2,:),...
                                       position(:,1),position(:,2) );

                    % create a new label for each new additional label
                    handles.labels(inpts) = maxID + 1;
                    maxID = maxID + 1;
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
        set( handles.refineClusts,'enable','on' );
        
        % update
        handles = update_sortmodel( handles );
        handles.clusterIDs = unique( handles.labels );
        handles.plotcolor = update_plot_colors( handles );
        st_update_plots( handles );
        guidata( hObject,handles );


    function measureQuality_Callback( hObject,eventdata,handles )
        % measures the quality of the clustering using the method specified by
        % handles.sortOptions.clusterMetric

        % check for projections
        if ~handles.availableProjection
            disp( 'No projections available' );
            return
        end

        % check for sorting
        if all( handles.labels == 0 )
            disp( 'Please sort clusters first' );
            return
        end

        % compute cluster quality and store into handles
        handles.R.clusterQuality...
         = measure_cluster_quality( handles.projection,...
                                    handles.labels,...
                                    handles.sortOptions.clusterMetric,...
                                    'euclidean'); % <-- to update, JS
        
        % plot the quality if the figure is available
        if ~isa( handles.qualityplot,'double' ) && isvalid( handles.qualityplot )
            st_plot_clusterQuality( handles );
        end
     
        guidata( hObject,handles );


    function handles = update_sortmodel( handles )
        % updates the sorting R for non-spectral sorting
        
        % check for model type
        switch handles.R.sortMethod
            case {'DBSCAN','Spectral'}
                % don't do anything for now, since model doesn't use probabilities
                return 

            case 'HDBSCAN'
                return
                
            otherwise
                % compute the means/covariances of the clusters
                proj = gather( handles.projection );
                model = get_sorting_model( proj,handles.labels,handles.R.sortMethod );
                
                % store the model and probabilities. Even if only a subset of
                % points are re-sorted, the probabilities for all points and
                % clusters will be recomputed. This provides a working model for
                % future spike sorting routines
                handles.R.sortModel = model;
                if isfield( model,'probabilities' )
                    handles.probabilities = model.probabilities;
                end
        end
                    
        
    %% random callbacks/fcns not assigned 
    function Untitled_1_Callback(hObject, eventdata, handles)