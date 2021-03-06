# mNEO
A MATLAB interpretation of [NEO](http://neuralensemble.org/neo/), a python package for organizing large datasets of multi-channel extracellular recordings.

mNEO is an all-purpose, modular package for analyzing and organizing large-scale electrophysiology recordings. While other excellent packages exist for a similar purpose (see "MClust", "KlustaKwik", "Spike2", and "SimpleClust" for some examples), each has its own shortcomming (including strict dependence on specific file formats, limited sorting routines, and lack of an integrated heirarchical organization system) that the mNeo package is aimed at resolving. This package is built completely in the MATLAB language, meaning the only pre-processing needed for interfacing with the package is loading one's data into the MATLAB workspace.

## Structure
The strength of this package comes from the combination of raw data preprocessing, spike detection and sorting, and an extensive multitiered system for organizing and storing data (highly inspired by the Neo Python package). The various classes used for data processing and organization are subclasses of the MATLAB *handle* class, meaning data and objects are all passed by reference. This results in a much more efficient data processing pipeline as objects (i.e. neurons, spikes, electrodes, etc.) and their corresponding data (i.e. spike waveforms, raw signals, etc.) grow in size. However, those unfamiliar with referenced data (as in Python) should familiarize themselves to avoid overwriting data accidentally. (Please also see the [heirarchy scheme](images/MatlabNeo_schematic.pdf))

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
The *Electrode* class is a container for *Signal* objects recoreded from the same, single channel. The *Electrode* object can contain multiple *Signal* objects (i.e. each pointing to a different *Epoch*), and can also belong to multiple *ChannelIndex* objects (as in a dense multi-electrode array where K-neighboring channels are considered a group). Note that while an *Electrode* object contains methods for pulling out the spike waveforms and raw voltage traces associated with that electrode, it cannot perform spike sorting or detection as, in general, such functionality is better suited to parallelization across multiple channels. Thus, one must ensure an *Electrode* is associated with 1 or more *ChannelIndex* objects for spike analysis.

*Parents:* ChannelIndex, Block

*Children:* Signal

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
Finally, the abstract *Container* class is not meant to be interfaced with directly, but provides the needed functionality for the classes listed above to reference one another via the "parent-child" framework. Because the *Container* class inherets from the builtin *handle* class, the above classes too inheret from the *handle* class, hence their ability to be reference instead of copied. All of the above classes can reference / de-reference one another through the `self.addChild`, `self.addParent`, `self.removeChild`, and `self.removeParent` methods, where `self` is replaced by the variable pointing to a particular class instance. Additionally, one may point local variables to particular children or parents of an object using the `self.getChild( 'child_class_type' )` and `self.getParent( 'parent_class_type' )`, where `child_class_type` and `parent_class_type` refer to the actual names of the classes contained within that particular instance. 

## Credit
While much of the package was written by me, certain features were either inspired by or taken directly from others. In particular:

1. Laurens van der Maaten's wonderful [dimension reduction toolbox](https://lvdmaaten.github.io/drtoolbox/)
2. The [fastICA](https://research.ics.aalto.fi/ica/fastica/) implementation of ICA
3. Joshua Stough's [select data function](http://www.mathworks.com/matlabcentral/fileexchange/37956-select-data) (used in the sortTool GUI)
4. Mo Chen's [EM-GMM and EM-VB](https://github.com/PRML/PRMLT) implementations
5. Yarpiz [DBSCAN](www.yarpiz.com)
6. Ulrike von Luxburg's [Spectral Clustering](https://github.com/Jorsorokin/neo-matlab/source/Functions/clustering/SpectralClustering.m)
7. Shimazaki's [kernel bandwidth estimation](https://github.com/Jorsorokin/neo-matlab/source/Functions/spike_analysis/sskernel.m)
8. Inspiration from the [Neo](http://neuralensemble.org/neo/) python package

## Uses 
For further information on using this package for data extraction and organization, see the [various examples](docs/examples), and the detailed [GUI documentation](docs/DRST_manual.pdf)

## To Do:
- [X] method for eliminating common spikes detected across multiple ChannelIndex objects
- [X] method for estimating best kernel for firing rates, and/or [adaptive kernels](https://link.springer.com/article/10.1007/s10827-009-0180-4)
- [X] functions for clustering via Variational Bayes and EM-TMM in the sortTool GUI
- [ ] method for visualizing the current heirarchy in the Block object (i.e. [connected graph](http://mathworld.wolfram.com/ConnectedGraph.html))
- [ ] dealing with different sorting parameters when updating the sortModel property of the Neuron object
- [ ] function for online sorting / assigning new spikes to current Neuron objects [wMSC.m](source/Functions/clustering/wMSC.m) & [HDBSCAN](https://github/Jorsorokin/HDBSCAN)
- [ ] interfacing with SQL for memory-efficient loading/saving and smart querying (SQLite?)
- [X] incorporating HDBSCAN into the sortTool GUI (see [HDSBCAN repo](https://github.com/Jorsorokin/HDBSCAN))
- [X] incorporating automatic cluster refinement & cluster quality measurement into sortTool GUI
- [ ] automatic rejection of clusters via cluster-quality thresholding (let users grab the threshold and slide it)
- [ ] padding of data outside of the range of the epochs with NaNs
- [ ] single-electrode referencing for spike detection 
- [ ] easily incorporated user-defined function for aligning action potentialss
