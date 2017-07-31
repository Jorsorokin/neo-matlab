# neo-matlab
A MATLAB-interpretation of [NEO](http://neuralensemble.org/neo/), a python package for organizing large datasets of multi-channel extracellular recordings.

Neo-matlab is an all-purpose, modular package for analyzing and organizing large-scale electrophysiology recordings. While other excellent packages exist for a similar purpose (see "MClust", "KlustaKwik", "Spike2", and "SimpleClust" for some examples), each has its own shortcomming (including strict dependence on specific file formats, limited sorting routines, and lack of an integrated heirarchical organization system) that the neo-matlab package is aimed at resolving. 

This package is built completely in the MATLAB language, meaning the only pre-processing needed for interfacing with the package is loading one's data into the MATLAB workspace. Although relying on the MATLAB language means this package has certain limitations (such as speed hits and more-rigid class structures compared to Python), MALTAB continues to be the defacto language for neural data analysis. Additionally, there are a plethora of thid-party functions developed by the MATLAB community that facilitate loading data from most systems (Plexon, Axona, NeuroNexus, OpenEphys...) and performing a wide range of sophisticated analysis. 

## Structure
The strength of this package comes from the combination of raw data preprocessing, spike detection and sorting, and an extensive multi-tiered system for organizing and storing data (highly inspired by the Neo Python package). The various classes used for data processing and organization are subclasses of the MATLAB *handle* class, meaning data and objects are all passed by reference rather than copy. This results in a much more efficient data processing pipeline as objects (i.e. neurons or channels) and their corresponding data (i.e. spikes or raw signals) grow in size. However, those unfamiliar with referenced data (as in Python) should familiarize themselves to avoid overwriting data accidentally. 

![alt text](https://github.com/Jorsorokin/neo-matlab/images/scheme.png "Heirarchy scheme")
**Figure 1: class heirarchy**

### Block 
At the top level, the *Block* class contains metadata (file names, paths, recording data, recording condition...) for a particular recording file, as well as *ChannelIndex* and *Epoch* classes as children. For the best organization, a data pipeline for each recording should begin with the creation of a *Block* instance, that accumulates the various channels, signals, and neurons discovered in the raw data. 

### ChannelIndex
A class for storing a group of channels together that logically should be contained as one object. For instance, the four channels of a [tetrode](https://en.wikipedia.org/wiki/Tetrode_(biology)) can be associated with one *ChannelIndex* instance, resulting in the subsequenct association of any action potential or neuron detected on any of these channels with the entire tetrode.

### Epoch
A sister class to the *ChannelIndex*, the *Epoch* class contains raw/processed signals and spikes associated with a particular window of time. An *Epoch* can be any user-defined window of time surrounding some event, such as the onset of a stimulus, the initiation of a decision go-cue, etc. Each such event can be stored into its own *Epoch* instance, which greatly facilitates trial-based analyses such as peri-stimulus time histograms (PSTHs), spike-triggered averaging (STA), or paired statistics between experimental events. 

### Signal  