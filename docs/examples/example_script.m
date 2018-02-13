% mNEO Example
% Jordan Sorokin, 2/11/2018
% For the Giocomo lab

%% Part 1: building a hierarchy from a multi-channel recording

% (a) create a block object to start the hierarchy
filePath = 'Z:\matlabscripts\SpikeSorting\example_data.mat';
name = 'mouse1_rec1';
date = '2017-02-20';
condition = 'baseline';
fs = 30000; % the sampling rate (fs = "frequency of sampling")
block = Block( name,date,condition,filePath );
block.print(); % print's some info to the screen about this block/recording

% (b) load in the data and create our Epochs and Electrodes. We will pretend that 
% we had 2 epochs, and have extracted +/- 5 seconds of data around each.
% Each column of the data matrix is one channel
epochTimes = [12, 28]; % seconds 
epochWindow = 5; % seconds, to extract around the epochTimes
load( filePath ); % load's in the data matrix
timestamps = (1:size( data,1 ))/fs;
inds = arrayfun( @(x)(timestamps > x-5 & timestamps <= x+5),epochTimes,'un',0 ); % indicies to extract for each epoch

% (c) create the Epochs, ChannelIndex, Electrodes, and Signals
epoch(1) = Epoch( epochTimes(1)-epochWindow,epochTimes(1)+epochWindow,1 ); 
epoch(1).name = 'epoch 1'; % you could assign any unique name that you wanted

epoch(2) = Epoch( epochTimes(2)-epochWindow,epochTimes(2)+epochWindow,2 );
epoch(2).name = 'epoch 2'; % ditto

% loop over channels
nChans = size( data,2 );
for ch = 1:nChans
    electrode(ch) = Electrode( ch ); % the electrode number
    electrode(ch).name = 'this is an electrode'; % optional name again
    
    
    for ep = 1:numel( epochTimes )
        
        % pull out data for epoch(ep) and data(ch)
        signal = Signal( data(inds{ep},ch),fs );
        
        % add this signal to the correspoinding epoch and electrode parent
        epoch(ep).addChild( signal );
        electrode(ch).addChild( signal );
    end
    
end

% now let's assume all the electrodes came from a silicon probe (in fact,
% they did), so we'll group them all into one single ChannelIndex parent
channelindex = ChannelIndex( 1 ); 
channelindex.addChild( electrode ); % note, we can add ALL of the electrodes at once. In fact, we can do this for any object
channelindex.name = 'shank 1'; % first shank of the silicon probe. Could also call this "tetrode 1" if using tetrodes

% finally add to the block and update the block. You will see that its
% parameters have changed, indicating the # of channelindex, electrodes,
% epochs, signals, and neurons it now has
block.addChild( epoch );
block.addChild( channelindex );
block.addChild( electrode );
block.print(); % calling block.print() also calls block.update()

% at this point, we could write our block to disk for safe keeping. To do
% so, you would provide an output directory, and call "block.write( outdir )".
% By default, the block oject writes to a file with the same name as
% "block.filename", and appended with "_extractedData.mat"

%% Part 2: detecting spikes and sorting using the GUI

% we will use the "double_flood_fill.m" algorithm, developed for
% high-density electrode arrays. It uses a flood-fill algorithm to find
% spikes specifically localized in BOTH time and space, while traditional
% detection is simply in time. For a tetrode, where each group of
% 4-channels is fairly well isolated from every other group, the
% traditional spike-detection may work equally well, and likely faster.

% (a) detect spikes and inspect the results
spikeThresh = 4; % x SD of the background
artifactThresh = 700; % uV...likely to be noise / movement
useMaskedDetection = true; % if false, then normal spike detection is used
channelindex.detectSpikes( spikeThresh,artifactThresh,useMaskedDetection );

% now let's use "block.print()" to see the results
block.print()

% if all worked well, your block object should now show 1 unit has been
% added. This "unit" represents an unsorted collection of Spike objects,
% one for each Electrode and Epoch. Let's take a look at it
neuron = channelindex.getChild( 'Neuron' );
disp( neuron );

% you can appreciate just how quickly spike counts increase with many
% channels. We've detected ~ 7000 spikes in just a 20 sec recording, using
% only a quarter of the total # of electrodes actually used to record this
% data set!

% (b) sort all of the spikes using the spike sorting GUI. One can sort both
% with and without the GUI, but visually inspecting sorting results often
% improves sorting quality.
projectionMethod = 'PCA';
nDimsToKeep = 10;               % keep the first 10 PCs
sortingMethod = 'HDBSCAN';      % hierarchical dbscan, download from github.com/jorsorokin/HDBSCAN
rejectThresh = 0.95;            % the closer to 1, the fewer the # of points rejected as noise (for HDBSCAN)
nNeighbors = 2;                 % a hyperparameter of HDBSCAN, check out the docs
minClustSize = 30;              % 30 spikes min to be considered a cluster
useGUI = true;

channelindex.sortSpikes( 'projMethod',projectionMethod,...
                         'level',nDimsToKeep,...
                         'sortMethod',sortingMethod,...
                         'nNeighbors',nNeighbors,...
                         'minclustsize',minClustSize,...
                         'useGUI',useGUI );
                     


% (c) resample the raw data since we don't need high sampling-rate as we've
% detected spikes already. This helps save space when saving the block
% object to disk!!!
p = 1;
q = 30; % resamples the signal by p/q...so sampling rate becomes 1 kHz
for ch = 1:block.nElectrodes
    block.getChild( 'Electrode',ch ).resampleSignals( p,q );
end

% (d) Print out the block again, and see how many neurons are now in the ChannelIndex object
block.print();

%% Part 3: visualization and simple analyses

% Now that we've isolated some units, let's look at some of their behaviors
% for the two epochs. Note that any neuron with ID = 0 means those spike
% waveforms that could not be sorted (i.e. multi-unit)
neurons = block.getNeurons(); % convenience method to get all Neurons in the block

% RASTERS
figure;
for j = 1:6
    subplot( 2,3,j )
    neurons(j+1).raster(5,1,2); % 1 secs before, 2 secs after eventOnset
end

% PSTH
neurons(4).psth(5,2,4,40);

% CUMULATIVE RASTER
figure; hold on
for j = 2:numel( neurons )
    [~,sptm] = neurons(j).getSpikes( 1 ); % pull out the spike times for first epoch only
    nSp = numel( sptm );
    scatter( sptm,ones( 1,nSp )*j,'.' ); % plot raster of spike times
end
ylabel( 'neurons' );
xlabel( 'time (s)' );
title( 'all neurons, epoch 1' );
darkPlot(gcf); % change the format of the figure
set( gca,'ylim',[1,block.nUnits+1] );

% RAW SPIKE WAVEFORMS
figure;
channelindex.plotSpikes( 1,[4,5] ); % first epoch, and neurons with IDs = 4,5 

% FEATURE MATRIX
neurons(8).plotFeatures();

% PLOT DATA FROM ELECTRODES
figure;
for j = 1:2
    subplot( 1,2,j );
    block.getChild( 'Epoch',j ).plotSignals();
end
