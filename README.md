# neo-matlab
A MATLAB interpretation of [NEO](http://neuralensemble.org/neo/), a python package for organizing large datasets of multi-channel extracellular recordings.

Neo-matlab is an all-purpose, modular package for analyzing and organizing large-scale electrophysiology recordings. While other excellent packages exist for a similar purpose (see "MClust", "KlustaKwik", "Spike2", and "SimpleClust" for some examples), each has its own shortcomming (including strict dependence on specific file formats, limited sorting routines, and lack of an integrated heirarchical organization system) that the neo-matlab package is aimed at resolving. 

This package is built completely in the MATLAB language, meaning the only pre-processing needed for interfacing with the package is loading one's data into the MATLAB workspace. Although relying on the MATLAB language means this package has certain limitations (such as speed hits and more-rigid class structures compared to Python), MALTAB continues to be the defacto language for neural data analysis. Additionally, there are a plethora of third-party functions developed by the MATLAB community that facilitate loading data from most systems (Plexon, Axona, NeuroNexus, OpenEphys...) and performing a wide range of sophisticated analysis. 

## Structure
The strength of this package comes from the combination of raw data preprocessing, spike detection and sorting, and an extensive multitiered system for organizing and storing data (highly inspired by the Neo Python package). The various classes used for data processing and organization are subclasses of the MATLAB *handle* class, meaning data and objects are all passed by reference. This results in a much more efficient data processing pipeline as objects (i.e. neurons, spikes, electrodes, etc.) and their corresponding data (i.e. spike waveforms, raw signals, etc.) grow in size. However, those unfamiliar with referenced data (as in Python) should familiarize themselves to avoid overwriting data accidentally. 

#### Block 
At the top level, the *Block* class contains metadata (file names, paths, recording data, recording condition...) for a particular recording file, and has methods for updating, printing, and writing all data contained within the block to disk. For the best organization, a data pipeline for each recording should begin with the creation of a *Block* instance that accumulates the various channels, signals, electrodes, and neurons discovered in the raw data. 

*Parents:* n/a

*Children:* ChannelIndex, Epoch, Electrode

#### ChannelIndex
A class for storing a group of channels together that should logically be contained as one object. For instance, the four channels of a [tetrode](https://en.wikipedia.org/wiki/Tetrode_(biology)) can be associated with one *ChannelIndex* instance, resulting in the subsequenct association of an action potential or neuron detected on any of these channels with all four channels of the tetrode. *ChannelIndex* objects contain methods for spike detection, sorting, and for interfacing with the *sortTool* GUI (see GUI folder). Although not all recording conditions may involve grouped electrodes, one should still associate spatially isolated electrodes with their own *ChannelIndex* instances for subsequent spike detection and sorting, as Electrode objects (below) do not have these capabilities.

*Parents:* Block

*Children:* Electrode, Neuron

#### Epoch
A sister class to the *ChannelIndex*, the *Epoch* class organizes raw or processed signals and spikes associated with a particular window of time. An *Epoch* can be any user-defined window of time surrounding some event, such as the onset of a stimulus, the initiation of a decision go-cue, etc. Each such event can be stored into its own *Epoch* instance, which greatly facilitates trial-based analyses such as peri-stimulus time histograms (PSTHs), spike-triggered averaging (STA), or paired statistics between experimental events. 

*Parents:* Block

*Children:* Signal, Spikes

#### Electrode
The *Electrode* class is a container for tying *Signal* and *Spikes* objects recoreded from the same, single channel. The *Electrode* object can contain multiple *Spikes* and *Signal* objects (i.e. each pointing to a different *Epoch*), and can also belong to multiple *ChannelIndex* objects (as in a dense multi-electrode array where K-neighboring channels are considered a group). Note that while an *Electrode* object contains methods for pulling out the spike waveforms and raw voltage traces contained in their *Spikes* and *Signal* children, it cannot perform spike sorting or detection as, in general, such functionality is better suited to parallelization across multiple channels. Thus, one must ensure an *Electrode* is associated with 1 or more *ChannelIndex* objects for spike analysis.

*Parents:* ChannelIndex, Block

*Children:* Signal, Spikes

#### Neuron
A *Neuron* object represents a (putatively) isolated neuron and its associated *Spikes* children (spike times and voltage waveforms). *Neuron* objects are assigned unique identifiers for each *ChannelIndex* independently. Thus, more than 1 *Neuron* with a certain ID may exist within the entire *Block* framework, however no two *Neuron* objects will have identical IDs and the same *ChannelIndex* parent. *Spikes* across all *Epochs* are associated with a particular *Neuron*, and under certain spike-sorting routines (i.e. EM-GMM, EM-TMM), that *Neuron* object will contain a working model of the low-dimensional spike-waveform distribution, allowing for the addition of future spikes.   

*Parents:* ChannelIndex

*Children:* Spikes

#### Signal
Raw and/or processed, continuous voltage traces are associated with a particular *Signal* object, itself contained within an *Electrode*. Multiple *Signal* objects within a single *Electrode* may exist, as in the case of multiple segments of data extracted from the entire recording. A *Signal* contains no children, and has functions including filtering, resampling, and noise estimation. By containing continuous data within a *Signal* object, one may reference or pass such data freely between variables without the overhead of copying large vectors.

*Parents:* Electrode, Epoch

*Children:* n/a

#### Spikes
Like the *Signal* class, the *Spikes* class has no children, but contains the raw voltage waveforms and times of action potentials detected from sister *Signal* objects. These *Spike* objects are contained within parent *Neuron* and *Epoch* objects, which facilitates clean organization and efficient retreival and analysis of neural firing patterns. *Spikes* are automatically created when running the spike detection method of the *ChannelIndex* object, provided there are *Signal* and *Electrode* objects to detect spikes from. Further, unsorted *Spike* objects will be automatically deleted and re-instantiated with the appropriate sorted waveforms belonging to new *Neuron* parents after running the spike-sorting method of the *ChannelIndex* object. Thus, one needn't worry about keeping track of which spike waveforms and times belong to particular neurons as the list of neurons grows, as this is taken care of intuitively and automatically.

*Parents:* Neuron, Epoch

*Children:* n/a

#### Container
Finally, the abstract *Container* class is not meant to be interfaced with directly, but provides the needed functionality for the classes listed above to reference one another via the "parent-child" framework. Because the *Container* class inherets from the builtin *handle* class, the above classes too inheret from the *handle* class, hence their ability to be reference instead of copied. All of the above classes can reference / de-reference one another through the "self.addChild", "self.addParent", "self.removeChild", and "self.removeParent" methods, where "self" is replaced by the variable pointing to a particular class instance. Additionally, one may point local variables to particular children or parents of an object using the "self.getChild( 'child_class_type' )" and "self.getParent( 'parent_class_type' )", where "child_class_type" and "parent_class_type" refer to the actual names of the classes contained within that particular instance. 

## Uses
Suppose we have recorded some neural activity using two tetrodes, during which two different stimuli were delivered. We wish to examine the effect of the stimulation on neural firing rates. First let us create an instance of the *Block* class to simplify subsequent data storage.

``` matlab
fileName = '07-21-17_mouse1_trial2.rec';
date = '07-21-17_15-45-13'; % date and time
condition = 'two stims';
filePath = '\Documents\Recordings\Tetrodes';

% create our Block instance. 
% For interpretability, it is easiest to assign the instances to variables with similar names as the classes they intantiate
block = Block( fileName,date,condition,filepath ); 
```

Now that we have created our Block, we want to next extract the relevant Signals. First lets load the raw data and create the two Epochs. 

``` matlab
% load the raw data
[data,samplingRate] = my_data_loader( [block.filepath,'\',block.filename] );
tetrodeChans = [1,2,7,8; 3,4,5,6]; % each row represents channels for one tetrode

% create our Epochs given the stimuli
stims = [20,43]; % 20s and 43s into the recording
stimWindow = 5; % extract 5s pre/post each stim
for i = 1:numel( stims )
    block.addChild( Epoch( stims(j)-stimWindow,stims(j)+postWindow,i ) ); % creates two Epoch instances, stored into the block parent
    thisEpoch = block.getChild( 'Epoch',i ); % point to the ith epoch
    thisEpoch.name = sprintf( 'stim %i' ); % gives it a specific name for us to remember later
    thisEpoch.addEvent( stims(i) ); % adds the time of the stim for each 
end
```

Note that we did not have to re-add the epochs back to the parent block after updating some of their properties. This is because the local variable "thisEpoch" pointed to the same physical memory as that particular epoch instance contained within the block. In other words, updating "thisEpoch" updated the underlying representation of memory pointed to by both "thisEpoch" and "block.getChild( 'Epoch',i )". This is the beauty of referencing vs. copying: we can take the viewpoint of the underlying memory through different objects, depending on what makes sense at the time. We also avoid having to re-copy all of the contents of the class instances after each update.

Now extract the actual voltage traces for each epoch and for each electrode. Then we will store back into the block object so that everything is updated.

``` matlab
% pull out the child epochs and get the starting and ending indices for pulling out raw voltage
epochs = block.getChild( 'Epoch' ); % pull out the epoch children
epochStarts = round( [epochs.startTime] * samplingRate );
epochEnds = round( [epochs.stopTime] * samplingRate );

% now extract the voltage traces for each Electrode (i.e. each column of the data) & epoch
for ch = 1:size( data,2 )
    electrode(ch) = Electrode( ch ); % create an Electrode object for this channel
    for ep = 1:numel( epochs )
        electrode(ch).addChild( Signal( data(epochStarts(ep,ch):epochEnds(ep,ch),samplingRate) ) );
        electrode(ch).getChild( 'Signal',ep ).filter( 300,6000 );
    end
end
block.addChild( electrode ); % adds the electrodes to the block
```

We just created 16 new Signal instances, each of which contains a *copy* of a small segment of data for a particular epoch and particular electrode. We then filtered the data of each signal between [300, 6000] Hz, and finally stored all the signals for a particular electrode into an Electrode instance. Note that our "electrode" variable is 1x8 here, where each dimension represents another electrode. Rather than performing "block.addChild( electrode )" for each pass of the for loop, we can just add all 1x8 electrode objects simultaneously. 

Finally, group our electrodes into two distinct tetrodes, so that we can perform spike detecting separately for the two groups of electrodes. We will detect spikes and leave them unsorted for now. This results in one Neuron instance created for each tetrode with ID = 0. In general, any Neuron with ID = 0 signifies unsorted spikes. When sorting, one may deem some spike waveforms as unsortable, in which case they will be added to the ID = 0 neuron. 

``` matlab
% create two ChannelIndex objects to group the electrodes into two tetrodes
for chGroup = 1:size( tetrodeChans,1 )
    channelindex = ChannelIndex( chGroups );
    channelindex.addChild( electrode(tetrodeChans(chGroup,:)) ); % adds appropriate Electrode children
    channelindex.name = sprintf( 'tetrode %i',chGroup ); % gives the channelindex a name for identifying which tetrode
    block.addChild( channelindex )
ende

% finally, detect spikes in our filtered data using the "channelindex.detectSpikes()" method
% each ChannelIndex object will now contain one Neuron instance, and that neuron will have two Spikes instances (one for each epoch)
block.getChild( 'ChannelIndex',1 ).detectSpikes( 5,1000 ); % spike threshold of 5 * noise SD, artifact threshold = 1000 uV
block.getChild( 'ChannelIndex',2 ).detectSpikes( 3,1000 ); % here we use a lower threshold, as this tetrode was less noisy
```

Now that we have created our electrodes & tetrode groups, extracted epochs and signals, and detected spikes, we want to save our block to disk for later analysis

``` matlab
saveDirectory = '\Documents\Analysis\Tetrodes';
block.write( saveDirectory );
block.print();
```

The "block.write()" method takes an output directory as an argument and saves the block object (and its children) under the filename: "block.filename_extractedData.mat". The "block.print()" method prints some useful summary information to the matlab screen. 








