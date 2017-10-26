function fig = st_create_optionsFig()
%
% creates the uifigure 'fig' for sorting options that is invoked 
% when clicking on the "options" menu of the main sorting tool figure

% create the tabgroup
set( 0,'DefaultFigureWindowStyle','normal' );
fig = figure( 'Visible','off','ToolBar','none','MenuBar','none',... 
	'Units','Normal','Position',[0.5 0.5 0.3 0.35],'Resize','off' );

% create individual uipannel
panelNames = {'General','Graphs'};
panel = zeros( 1,numel(panelNames ) );
positions = {[0.02 0.02 0.4 0.98],[0.45 0.02 0.53 0.98]};
for i = 1:numel( panelNames )
	panel(i) = uipanel( 'Parent',fig,'Title',panelNames{i},...
		'Tag',panelNames{i},'Position',positions{i} );
end

%% General panel
names = {'Search for cluster #','Rejection probability:','0.05','Outlier threshold:','0.9',...
    'Measure cluster dispersion','Default # of clusters:','2' };
types = {'checkbox','text','edit','text','edit','checkbox','text','edit'};
tags = {'searchForK','rejectProb','probBox','outlierThresh','outlierBox','clustErr','clusterNum','kBox'};
positions = {[0.05 0.85 0.7 0.06],...                       % searchForK
             [0.05 0.65 0.5 0.06],[0.57 0.65 0.15 0.06],... % rejectProb + probBox
             [0.05 0.45 0.5 0.06],[0.57 0.45 0.15 0.06],... % outlierThresh + outlierBox
             [0.05 0.75 0.7 0.06],...                       % clustErr
			 [0.05 0.55 0.5 0.06],[0.57 0.55 0.15 0.06]};   % clusterNum + kBox

for j = 1:numel( names )
	uicontrol( 'Parent',panel(1),'Style',types{j},'String',names{j},...
			'Units','Normal','Position',positions{j},'Tag',tags{j} );
end

% align handles
h = findobj( panel(1),'Tag','rejectProb','-or','Tag','probBox' );
align( h,'VerticalAlignment','Center' );
clear h

%% Graph panel
names = {'kNN','mkNN','eps',...
	'# neighbors:','15','epsilon radius [0:1]:','0.1',...
    'RBF sigma:','1','minclustsize','5',...
    'Laplacian normalization','Kernel weighting'};
tags = {'','','',...
	'nNeighbors','nnBox','epsilon','epsBox','sigmaRBF','sigmaBox',...
	'minclust','minclustBox','laplaceNorm','weightGraph'};
types = {'Radio','Radio','Radio',...
	'text','edit','text','edit','text','edit',...
    'text','edit','checkbox','checkbox'};
positions = {[0.05 0.8 0.8 0.12],[0.05 0.5 0.8 0.12],[0.05 0.2 0.8 0.12],...% radio bttns
		[0.05 0.45 0.5 0.06],[0.6 0.45 0.15 0.06],...                       % nn box
        [0.05 0.35 0.5 0.06],[0.6 0.35 0.15 0.06],...                       % epsilon box
        [0.05 0.25 0.5 0.06],[0.6 0.25 0.15 0.06],...                       % sigma box
        [0.05 0.15 0.5 0.06],[0.6 0.15 0.15 0.06],...                       % minclust box 
		[0.48 0.8 0.5 0.08],...                                             % laplacian checkbox
        [0.48 0.7 0.5 0.08]};                                               % weightGraph box
defaults = {1,0,0,nan,10,nan,0.1,nan,1,nan,5,1,0};

% first create the radiobuttongroup
bg = uibuttongroup( panel(2),'Title','Neighborhood','Units','Normal',...
			'Position',[0.05 0.55 0.3 0.4],'Tag','neighborType' );

for j = 1:numel( names )
	if strcmp( types{j},'Radio' )
		uicontrol( 'Parent',bg,'Style',types{j},'String',names{j},...
			'Units','Normal','Position',positions{j},'Value',defaults{j},'Tag',tags{j} );
	else
		uicontrol( 'Parent',panel(2),'Style',types{j},'String',names{j},...
			'Units','Normal','Position',positions{j},'Value',defaults{j},'Tag',tags{j} );
	end
end

% align handles
h = findobj( panel(2),'Tag','nNeighbors','-or','Tag','nnBox' );
align( h,'VerticalAlignment','Center' );
h = findobj( panel(2),'Tag','epsilon','-or','Tag','epsBox' );
align( h,'VerticalAlignment','Center' );
h = findobj( panel(2),'Tag','sigmaRBF','-or','Tag','sigmaBox' );
align( h,'VerticalAlignment','Center' );
h = findobj( panel(2),'Tag','minclust','-or','Tag','minclustBox' );
align( h,'VerticalAlignment','Center' );
h = findobj( panel(2),'Tag','nNeighbors','-or','Tag','epsilon','-or','Tag','sigmaRBF','-or','Tag','minclust' );
align( h,'HorizontalAlignment','Right' );

end