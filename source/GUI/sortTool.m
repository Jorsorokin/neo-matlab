function varargout = sortTool(varargin)
    % SORTTOOL MATLAB code for sortTool.fig
    %
    % GUI for facilitating spike sorting of large databases. Inputs are optional,
    % as one can load in the various files through the GUI itself.
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
    %              observation, rows are variables)
    %
    %   ( times ) - the times of the spike waveforms (as one long
    %               vector, in seconds)
    %
    %   ( projection ) - the projections of the data onto lower dims
    %
    %   ( labels ) - the labels from a previous spike-sorting routine
    %
    %   ( trials ) - trial-labels of the spikes (i.e. one can upload a
    %                a large data matrix with spikes from different trials,
    %                but to get an accurate representation of ISI one needs
    %                information regarding which spikes are from which
    %                trial. Default is to assume the same trail).
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
    %                   W - the projection matrix
    %                   projMethod - the type of projection specified
    %                   sortMethod - the sorting method (kmeans, EM,...)
    %                   sortParams - covariances (sigma), means (mu), and
    %                                weights (w) of the sorting
    %                   probability - the probabilities of each spike
    %                                 belonging to each cluster
    %
    % Last edited by Jordan Sorokin, 7/23/2017

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

        % set up user-defined variables in the handles structure
        handles.figure1.Visible = 'off';
        handles.figure1.Name = 'Dimensionality Reduction & Sorting Tool';
        handles.data = nan;
        handles.times = nan;
        handles.projection = nan;
        handles.labels = nan;   
        handles.trials = nan;
        handles.selectedPoints = nan;
        handles.manualClust = {};
        handles.plotDims = [1,2];                       % default plotting 1st & 2nd columns
        handles.nDim = 3;                               % defaults to 3D projection 
        handles.model = struct; 
        handles.model.projectMethod = 'PCA';            % defaults to PCA first
        handles.model.sortMethod = 'EM-GMM';            % defaults to gaussian mixture modeling
        handles.model.W = nan;        
        handles.model.probabilities = nan;
        handles.model.keptPts = nan;
        handles.allPlotColor = [ [.8, .8, .8];...       % light gray        (0)
                                 [0.1, 0.74, 0.95];...  % deep sky-blue     (1)
                                 [0.95, 0.88, 0.05];... % gold/yellow       (2)
                                 [0.80, 0.05, 0.78];... % magenta           (3)
                                 [0.2, 0.7, 0.20];...   % lime green        (4)
                                 [0.95, 0.1, 0.1];...   % crimson red       (5)   
                                 [0.64, 0.18, 0.93];... % blue-violet       (6)
                                 [0.88, 0.56, 0];...    % orange            (7)
                                 [0.4, 1.0, 0.7];...    % aquamarine        (8)
                                 [0.95, 0.9, 0.7];...   % salmon-yellow     (9)
                                 [0, 0.4, 1];...        % blue              (10)
                                 [1, 0.41, 0.7];...     % hot pink          (11)
                                 [0.6, 1, 0];...        % chartreuse        (12)
                                 [0.6, 0.39, 0.8];...   % amtheyist         (13)
                                 [0.82, 0.36, 0.36,];...% indian red        (14)
                                 [0.53, 0.8, 0.98];...  % light sky blue    (15)
                                 [0, 0.85, 0.1]...      % forest green      (16)
                                ];
        % repeat colors in case we have many identified neurons
        handles.allPlotColor = [handles.allPlotColor; repmat( handles.allPlotColor(2:end,:),4,1 )]; 
        handles.plotcolor = nan;        
        
        % GET THE OPTIONAL INPUTS
        p = st_check_inputs( varargin );
        handles.data = p.data;
        handles.labels = p.labels;
        handles.times = p.times;
        handles.trials = p.trials;
        handles.projection = p.projection;
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
        if ~isnan( handles.projection )
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
        if ~isa( handles.waveformplot,'double' ) && isvalid( handles.waveformplot )
            close( handles.waveformplot.Parent );
        end
        if ~isa( handles.isiplot,'double' ) && isvalid( handles.isiplot )
            close( handles.isiplot.Parent );
        end
        if ~isa( handles.rasterplot,'double' ) && isvalid( handles.rasterplot )
            close( handles.rasterplot.Parent );
        end
        
                 
    % --- clears the current figure
    function CloseMenuItem_Callback(hObject, eventdata, handles)
        % hObject    handle to CloseMenuItem (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                             ['Close ' get(handles.figure1,'Name') '...'],...
                             'Yes','No','Yes');
        if strcmp(selection,'No')
            return;
        end
        delete( handles.figure1 ); 
        close all
        

    % --- Executes on button press in savebutton.
    function savebutton_Callback(hObject, eventdata, handles)
        % closes the GUI   
        handles.figure1.Visible = 'off';
        
        
    % --- Executes on button press in delete
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

        
% -------------------------------------------------------------------------  
% loading & projections
% ------------------------------------------------------------------------- 
    % --- Executes on button press in loaddata.
    function loaddata_Callback(hObject, eventdata, handles)
        % hObject    handle to loaddata (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        st_clear_plots( handles )

        % allow user to load data from .mat, .csv, or .txt file
        [file,filepath] = uigetfile( {'*.mat';'*.csv';'*.txt'},...
            'select the data matrix to sort' );
        handles.data = load( [filepath,filesep,file] );
        handles.labels = zeros( 1,size( data,2 ),'uint8' );
        handles.plotcolors = update_plot_colors( handles );
        handles.selectedPoints = false( 1,numel( handles.labels ) );
        handles.model.keptPts = single( 1:size( handles.data,2 ) );

        % update the GUI
        guidata( hObject, handles );
        
        
    % --- Executes during object creation, after setting all properties.
    function mappingmethod_CreateFcn(hObject, eventdata, handles)
        % hObject    handle to mappingmethod (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    empty - handles not created until after all CreateFcns called
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
             set(hObject,'BackgroundColor','white');
        end

        % set the available options
        set( hObject, 'String',...
            {'PCA (principal components analysis)',...
            'kPCA (kernel PCA)',...
            'pPCA (probabilistic PCA',...
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
        % hObject    handle to mappingmethod (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        allmethods = {'PCA','KernelPCA','ProbPCA','ICA','LLE','HessianLLE','NPE','SPE',...
                         'MVU','FastMVU','Laplacian','DiffusionMaps','SNE','tSNE',...
                         'Isomap','Sammon','Autoencoder'};
        value = get( hObject,'value' );
        handles.model.projectMethod = allmethods{value};

        % update the GUI
        guidata( hObject, handles );

        
    % --- Executes on button press in projectdata.
    function projectdata_Callback(hObject,eventdata,handles)
        % hObject    handle to projectdata (see GCBO)
        % handles    structure with handles and user data (see GUIDATA)

        % project the data using the method in "handles.projectMethod"
        [handles.projection,handles.model.W] = st_project_data( handles );
        
        % now plot onto the main figure
        handles.plotDims( handles.plotDims > handles.nDim ) = handles.nDim;
        st_plot_projections( handles )
        handles.selectedPointPlot = st_plotSelectedData( handles );

        % update the GUI
        guidata( hObject, handles );

                
    % --- Executes during object creation, after setting all properties.
    function dim1_CreateFcn(hObject, eventdata, handles)
        % hObject    handle to dim1 (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    empty - handles not created until after all CreateFcns called

        % Hint: edit controls usually have a white background on Windows.
        %       See ISPC and COMPUTER.
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end
        set(hObject,'String',1); % automatically sets to 1 when initializing
        guidata( hObject,handles );

        
    function dim1_Callback(hObject, eventdata, handles)
        % hObject    handle to dim1 (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)

        % Hints: get(hObject,'String') returns contents of dim1 as text
        %        str2double(get(hObject,'String')) returns contents of dim1 as a double
        handles.plotDims(1) = str2double( get( hObject,'String' ) );
        guidata( hObject,handles );

        
    % --- Executes during object creation, after setting all properties.
    function dim2_CreateFcn(hObject, eventdata, handles)
        % hObject    handle to dim2 (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    empty - handles not created until after all CreateFcns called

        % Hint: edit controls usually have a white background on Windows.
        %       See ISPC and COMPUTER.
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end

        set(hObject,'String',2); % automatically sets to 2 when initializing
        guidata( hObject, handles );
        
        
    function dim2_Callback(hObject, eventdata, handles)
        % hObject    handle to dim2 (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)

        % Hints: get(hObject,'String') returns contents of dim2 as text
        %        str2double(get(hObject,'String')) returns contents of dim2 as a double
        handles.plotDims(2) = str2double( get( hObject,'String' ) );
        guidata( hObject,handles );     
        
               
    function nDim_CreateFcn(hObject, eventdata, handles)
        % hObject    handle to nDim (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    empty - handles not created until after all CreateFcns called

        % Hint: edit controls usually have a white background on Windows.
        %       See ISPC and COMPUTER.
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end
        set( hObject,'String',3 ); % default = 3
        guidata( hObject,handles );
        

    function nDim_Callback(hObject, eventdata, handles)
        % hObject    handle to nDim (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        handles.nDim = str2double( get( hObject,'String' ) );
        guidata( hObject,handles );    
  
        
% ------------------------------------------------------------------------- 
% misc. plotting procedures (updating colors, replotting, etc.)
% ------------------------------------------------------------------------- 
     % --- Executes on button press in replogt
    function replogt_Callback(hObject, eventdata, handles)
        % update the plots according to new labels or projection dims
        st_plot_projections( handles );
        handles.selectedPointPlot = st_plotSelectedData( handles );
        st_update_plots( handles );
        guidata( hObject,handles );
        
        
    % --- Executes on button press in selectdata.
    function selectdata_Callback(hObject, eventdata, handles)
        % allows user to select individual data points in the main plot     
        if isnan( handles.data )
            return;
        end

        % get current axes, clear previous selected points
        %handles = st_clearSelectedData( handles );

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
        if isnan( handles.data )
            return;
        end
        
        % get current axes, clear previous points
        %handles = st_clearSelectedData( handles );
        
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
        handles = st_clearSelectedData( handles );
        
        % update
        guidata( hObject,handles );
        
                    
    % --- updates plotting colors of all points
    function colors = update_plot_colors( handles )
        
        % check if any labels (i.e. any data)
        if isnan( handles.labels )
            colors = nan;
            return;
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
        if isa( handles.waveformplot,'double' ) || ~isvalid( handles.waveformplot )
            handles.waveformplot = st_create_waveplot( handles ); 
        end
            
        % check if visible; if not, plot
        switch handles.waveformplot.Parent.Visible
            case 'on'
                cla( handles.waveformplot );
                handles.waveformplot.Parent.Visible = 'off'; % hide the plot
            otherwise
                st_plot_waveforms( handles );
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
                
        
% ------------------------------------------------------------------------- 
% auto/manual sorting
% ------------------------------------------------------------------------- 
    % --- Initiates the sort-method dropdown list
    function sortmethod_CreateFcn(hObject,eventdata,handles)
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
             set(hObject,'BackgroundColor','white');
        end

        % set the available options
        set( hObject, 'String',...
            {'EM-GMM (gaussian mixture model)',...
            'EM-TMM (t-dist. mixture model)',...
            'Km (K-means)',...
            'VB (variational bayes)'} );

        % update the GUI
        guidata( hObject, handles ); 
        
        
    % --- Executes on button press in sortmethod
    function sortmethod_Callback(hObject,eventdata,handles)
        allmethods = {'EM-GMM','EM-TMM','Km','VB'};
        value = get( hObject,'value' );
        handles.model.sortMethod = allmethods{value};

        % update the GUI
        guidata( hObject, handles );
    
        
    % --- Executes on button press in autosort.
    function autosort_Callback(hObject, eventdata, handles)        
        warning off;
        
        % clear any data selection
        handles = st_clearSelectedData( handles );
        
        % first project the data 
        if isempty( handles.data )
            errordlg( 'No data available. Please load data!' );
            return;
        elseif isnan( handles.projection )
            disp( 'please project your data first!' );
            return;
        end
        
        % now sort via automatic cluster search or manual cluster number
        searchForK = questdlg( 'Search for number of clusters?','yes','no' );
        switch searchForK
            case 'Yes'
                [~,K] = findClustNum( handles.projection',2,15 );
            case 'No'
                K = ( inputdlg( 'Number of clusters' ) );
                K = str2double( cell2mat( K ) );
            case 'Cancel'
                return;
        end
        
        % now do the actual clustering using the method specified
        [labels,handles.model.sortModel,handles.model.probabilities] = ...
            st_sort_clusters( handles.projection,K,handles.model.sortMethod );
        
        % change the labels to contain the original label plus consecutive
        % labels for each new cluster
        labels = labels + min( handles.labels ) - 1; % subtract 1 to start with the minimum of the previous labels
        
        % refine the clusters
        switch handles.model.sortMethod
            case {'EM-GMM','EM-TMM','VB'}
                [handles.labels,keep]...
                    = st_refine_cluster( labels,handles.model.probabilities,0.1 );
                %handles = st_rmdata( handles,~keep ); 
        end
        
        % update 
        [handles.labels,handles.probabilities] = update_labels( handles );
        handles.plotcolor = update_plot_colors( handles );
        st_update_plots( handles );
        guidata( hObject,handles );
        
        warning on;
    
        
    % --- Executes on button press in manualsort.
    function manualsort_Callback(hObject, eventdata, handles)
        % hObject    handle to manualsort (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        
        % clear any data selection
        handles = st_clearSelectedData( handles );
        
        % check if data has been plotted
        axes( handles.mapplot );
        if isempty( handles.mapplot.Children )
            if ~isnan( handles.projection )
                st_plot_projections( handles );
            end
        end        
        
        % set visibility of the new buttons to "on"
        set( handles.newLabel,'enable','on','visible','on' );
        set( handles.mergeLabel,'enable','on','visible','on' );
        set( handles.deleteLabel,'enable','on','visible','on' );
        set( handles.finishManual,'enable','on','visible','on' );
        set( handles.autosort,'enable','off' );
        set( handles.projectdata,'enable','off' );
        set( handles.loaddata,'enable','off' );
        set( handles.replogt,'enable','off' );
        
        % update
        guidata( hObject,handles );
        
        
    % --- Executes on button press in newLabel.
    function newLabel_Callback(hObject, eventdata, handles)
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
        st_clearSelectedData( handles );
        guidata( hObject,handles );

        
    % --- Executes on button press in deleteLabel.
    function deleteLabel_Callback(hObject, eventdata, handles)
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

                    % TO DO
                    % make probability vector even without having first
                    % called the AUTOSORT function...
                    % create a new label for each new additional label
                    handles.labels(inpts) = max( handles.labels )+1;
                    if ~any( isnan( handles.model.probabilities ) )
                        handles.model.probabilities(:,end+1) = 0;
                        handles.model.probabilities(inpts,:) = 0;
                        handles.model.probabilities(inpts,end) = 1;
                    end
                end
                delete( handles.manualClust{i} );
            end
        end
        
        % now hide the extra buttons    
        set( handles.newLabel,'enable','off','visible','off' );
        set( handles.mergeLabel,'enable','off','visible','off' );
        set( handles.deleteLabel,'enable','off','visible','off' );
        set( handles.autosort,'enable','on' );
        set( handles.finishManual,'enable','off','visible','off' );
        set( handles.projectdata,'enable','on' );
        set( handles.loaddata,'enable','on' );
        set( handles.replogt,'enable','on' );
        
        % update
        [handles.labels,handles.probabilities] = update_labels( handles );
        [handles.model.sortModel.mu,handles.model.sortModel.sigma]...
            = get_cluster_description( handles.projection,handles.labels );
        handles.plotcolor = update_plot_colors( handles );
        st_update_plots( handles );
        guidata( hObject,handles );
    
        
    % --- Update all the labels
    function [labels,probabilities] = update_labels( handles )
        % updates the labels so that it adds new labels for each new clust
        labels = handles.labels;
        probabilities = handles.model.probabilities;
        uID = unique( labels(labels~=0) );
        possibleLabels = 1:max( uID );
        noProb = any( isnan( probabilities ) );
%         for i = 1:numel( uID )
%             prevLabels = ismember( labels,uID(i) );
%             if ~noProb                
%                 if any( prevLabels )
%                     probabilities(prevLabels,i) = probabilities( prevLabels,uID(i) );
%                 end
%             end
%             labels(prevLabels) = i;
%         end
        
        % remove extra probabilities not associated with any new label
        if ~noProb & max( labels ) < size( probabilities,2 )
            probabilities = probabilities(:,1:max( labels ));
        end
                    
        
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