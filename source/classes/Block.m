classdef Block < Container
    
    properties
        filename
        date
        condition
        filepath
        nEpochs = 0;
        nChanInds = 0;
        nUnits = 0;
        nSignals = 0;
    end
    
    methods
        
        function self = Block( filename,date,condition,filepat h)
            % self = Block( filename,date,condition,filepath )
            %
            % Initiate an instance of the Block class. A Block is 
            % the top-level container for a given recording session,
            % and contains ChannelIndex and Epoch classes as children.
            % 
            % The main purpose of the Block container is to tie the 
            % data (including channels, epochs, signals, and neurons)
            % together into a unified recording session.
            %
            % Children:
            %   Epoch
            %   ChannelIndex
            %
            % Parents:
            %   none
            %
            % Methods:
            %   print
            %   write
            %   
            %       * see also methods in the Container class

            self.filename = filename;
            self.date = date;
            self.condition = condition;
            self.filepath = filepath;
        end
        
        
        function print( self )
            % print( self )
            % 
            % displays the file metadata and children
            fprintf( '%s (%s)\nRecorded on %s\n',...
                self.filename,self.condition,self.recTime,self.date );
            ChInd = self.getChild( 'ChannelIndex' );
            Ep = self.getChild( 'Epoch' );
            nNeurons = 0;
            for c = 1:numel( ChInd )
                nNeurons = nNeurons + ChInd(c).nUnits;
            end
            fprintf( '%i identified units over %i ChannelIndex & %i Epoch children\n',...
                nNeurons, numel( ChInd ), numel( Ep ) );
        end
        
        
        function addChild(self,child)
            switch class( child )
                case {'Epoch'}
                    addChild@Container( self,child );
                    self.nEpochs = self.nEpochs + numel( child );
                case {'ChannelIndex'}
                    addChild@Container( self,child );
                    self.nChanInds = self.nChanInds + numel( child );
                otherwise 
                    error( 'Only ChannelIndex and Epoch objects are valid children' );
            end
        end
        

        function addParent(~,~)
            disp( 'Block objects have no parents' );
        end
        
        
        function write( self,outpath )
            % write( self,path )
            %
            % saves the Block and its children into the folder specified by
            % "outpath". The file will be saved as "self.filename.mat".
            
            if ~isdir( outpath )
                mkdir( outpath )
                addpath( outpath );
            end

            disp( 'Writing the Block to disk...' )
            save( [outpath,filesep,self.filename,'.mat'],self );
        end
                
    end
    
end

