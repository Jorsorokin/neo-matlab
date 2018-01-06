classdef Electrode < Container
    % self = Electrode( electrodeNum )
    %
    % Creates an instance of the Electrode class.
    % An Electrode is a container for multiple Signal objects that
    % originate from the same recording channel.
    %
    % An Electrode object can belong to multiple ChannelIndex
    % objects (as in the case of dense, multi-electrode 
    % probes/arrays), and can also have multiple Signal children,
    % which refer to different Epochs. 
    %
    % Children:
    %   Signal
    %
    % Parents:
    %   ChannelIndex
    %   Block
    %
    % Properties:
    %   nSignals - the number of channels contained in this Signal
    %   chanInd - the parent ChannelIndex numbers
    %   nChanInd - the number of ChannelIndex parents
    %   electrodeNum - a unique number refering to this electrode
    %   name - a user-defined tag for the electrode (i.e. "shank 3, chan 2")
    %   
    % Methods:
    %   plot
    %   getSpikes
    %   getVoltage
    %   resampleSignals
    %   
    %       * see also methods in the Container class  

    properties
        nSignals = 0;
        nChanInd = 0;
        chanInd = [];
        electrodeNum
        name
    end
    
    methods
        
        function self = Electrode( electrodeNum )            
            self.electrodeNum = electrodeNum;
        end
        
        
        function addChild( self,child )
            switch class( child )
                case 'Signal'
                    addChild@Container( self,child );
                    self.nSignals = self.nSignals + numel( child );
                    for j = 1:numel( child )
                        child(j).electrode = self.electrodeNum;
                        child(j).chanInd = [child(j).chanInd,self.chanInd]; % adds the associated ChannelIndex
                    end
                otherwise
                    error( 'Only Signal objects are valid children' );
            end
        end 


        function addParent( self,parent )
            switch class( parent )
                case 'ChannelIndex'
                    addParent@Container( self,parent );
                    parent.nSignals = parent.nSignals + self.nSignals;
                    parent.nChanInd = parent.nChanInd + 1;
                    self.chanInd(end+1) = parent.chanIndNum;
                case 'Block'
                    parent.nElectrodes = parent.nElectrodes + 1;
                otherwise
                    error( 'Only ChannelIndex objects are valid parents' );
            end
        end
        
        
        function [snips,times,epNum,neuronID,chIndNum] = getSpikes( self,varargin )
            % [snips,times,epNum,neuronID,chIndNum] = getSpikes( self, (epochs) )
            %
            % pull out the Spikes objects (if any) associated with
            % this Electrode. Can optionally specify which epochs to extract
            % spikes from. 
            %
            % Because an Electrode can belong to multiple ChannelIndex
            % objects, the "chaninds" output also signifies which spike-waveforms
            % and times come from which ChannelIndex parents. This can aid in removing
            % redundant spikes that were detected multiple times from the same 
            % electrode, but in different channel groups. 

            % check inputs
            if nargin > 1 && ~isempty( varargin{1} )
                epochs = varargin{1};
            else
                epochs = [];
            end

            snips = [];
            times = [];
            epNum = [];
            chIndNum = [];
            neuronID = [];
            
            % pull out neurons associated with this electrode (same channelindex)
            neurons = self.getSibling( 'Neuron','ChannelIndex' );
            if isempty( neurons )
                disp( 'No spikes available' );
                return
            end

            % loop over neuron objects, pull out spikes, times, epoch nums
            for i = 1:numel( neurons )
                [v,t,ep] = neurons(i).getSpikes( epochs );
                t = reshape( t,1,numel( t ) ); % makes into one long vector
                t(isnan(t)) = [];

                % get the chanind vector
                ch = repmat( neurons(i).chanInd,1,numel( ep ) );
                id = repmat( neurons(i).ID,1,numel( ep ) );

                % add to the outputs
                thisElectrode = ismember( neurons(i).getParent( 'ChannelIndex' ).chanIDs,self.electrodeNum );
                snips = [snips, squeeze( v(:,:,thisElectrode) )];
                times = [times,t];
                epNum = [epNum,ep];
                chIndNum = [chIndNum,ch];
                neuronID = [neuronID,id];
            end
        end


        function voltage = getVoltage( self,varargin )
            % voltage = getVoltage( self,(epochs) )
            % 
            % pull out the voltage waveforms from the Signal
            % children (if any) associated with this Electrode.
            % 
            % voltage is an n x m matrix, with n = max number of points 
            % across epochs, and m = # of signals (i.e. # of epochs).

            % check inputs
            if nargin > 1 && ~isempty( varargin{1} )
                epochs = varargin{1};
            else
                epochs = 1:self.nSignals;
            end

            % pull out the signal objects
            signals = self.getChild( 'Signal',epochs );
            if ~isempty( signals )

                % preallocate the voltage matrix
                n = max( [signals.nPoints] );
                voltage = nan( n,self.nSignals );

                % loop over signals
                for sig = 1:self.nSignals
                    voltage(1:signals(sig).nPoints,sig) = signals(sig).voltage;
                end
            else
                disp( 'No signals available' );
                voltage = nan;
            end
        end
        
        
        function resampleSignals( self,p,q )
            % resampleSignals( self,p,q )
            %
            % extract voltage traces from child Signal objects and resample
            % according to the ratio p/q
            for j = 1:self.nSignals
                self.getChild( 'Signal',j ).resample( p,q );
            end
        end
                
        
        function hL = plot( self )
            % hL = plot( self )
            %
            % plot signals in the current Electrode object and returns the handle.
            %
            % get the epoch start/end time if this Electrode object has an
            % Epoch partner. Plot relative to this start/stop time. Also
            % plot "." for spike times (if any), color coded by the Neuron ID
            if self.nSignals == 0
                disp( 'No signals available' );
                return
            end
            
            % get the epochs and create a time-matrix for plotting
            signals = self.getChild( 'Signal' );
            ep = self.getPartner( 'Epoch','Signal' );   
            if ~isempty( ep )
                events = [ep.eventTime];
            else
                events = [];
            end
            
            % get the voltage traces
            voltage = self.getVoltage();
            fs = signals(1).fs;

            % plot the signal on the full time-scale
            gca; hold on;
            [~,hL] = multisignalplot( voltage,fs );
            for i = 1:self.nSignals

                % plot the epoch event as a small red dot
                if ~isempty( events )
                    scatter( events(i)-ep(i).startTime,max( hL(i).YData ),80,'r.' );
                end
            end

            % clean up graph
            xlabel( 'time (s)' );
            ylabel( 'epoch' );
            title( sprintf( 'Electrode: %i, all Epochs',self.electrodeNum ) );  
            darkPlot( gcf );           
            hold off;
        end
        
    end % methods
    
end