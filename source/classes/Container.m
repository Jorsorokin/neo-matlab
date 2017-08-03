classdef Container < handle
    % Container < handle
    %
    % An abstract superclass for the classes: 
    % Block, Epoch, ChannelIndex, Neuron, Electrode, Signal, and Spikes 
    %
    % Provides the functionality for each of these objects to
    % reference one another via the "child-parent" syntax. A Container
    % should not be directly instantiated, but instead is automatically
    % inherited by the other objects listed above. 
    %
    % Because a Container inherets from "handle", and each of the above 
    % objects inherets from Container, everything is accessed by reference 
    % rather than copying. This can dramatically improve speed as memory 
    % isn't allocated when assigning a reference to one of these objects 
    % to a local variable, HOWEVER the programmer should be aware that 
    % references point to the same memory, and thus changing the value of 
    % a property of one reference to an object will change the value of 
    % that property in another reference to the same object.
    %
    % Children:
    %   none
    %
    % Parents:
    %   none
    %
    % Properties:
    %   parent - an empty (to begin) cell array. As parents are added, each
    %            new type will be stored as its own n-dimensional element,
    %            where the dimension "n" refers to the number of parents 
    %            of that type.
    %   children - ditto, but for added children.
    %   
    %       * in general, a "parent" is a higher-level object containing
    %         a lower-level child as a reference. Thus, "Neuron" is the
    %         parent of "Spikes", since each neuron should have its
    %         associated series of spike trains. 
    %
    % Methods
    %   addChild
    %   addParent
    %   removeChild
    %   removeParent
    %   getChild
    %   getParent
    %   getPartner
    %   getSibling
    %   deleteSelf <-- note, calling this will clear the current object from
    %                  memmory, as well as any reference to the object
    %
    %       * methods are intentionally redundant. I.e. if both a Neuron
    %       and its associated Spikes objects exist, then calling
    %       "Neuron.removeChild( 'Spikes' )" or "Spikes.removeParent(
    %       "Neuron" ) will both accomplish the same thing: the parent
    %       Neuron will de-reference the Spikes child, AND the Spikes child
    %       will de-reference the Neuron parent. The redundancy in calls
    %       gives this Container abstract class a highly modular form,
    %       allowing the programmer to take the viewpoint of whichever
    %       object (the child or the parent) seems more intuitive for
    %       his/her purposes at any given time.
    
    properties
        parent = {};
        children = {};
    end
    
    methods
        
        function addChild( self,child )
            % addChild( self,child )
            %
            % update the children. If this type already exists as 
            % a child, add to it by appending to the end
            
            % check if an empty version of this type of child
            childClass = class( child );
            [previous,ind] = self.getChild( childClass );
            
            if ~isempty( previous ) && isa( previous,childClass )
                self.children{ind}(end+1) = child;
            elseif isempty( previous ) && isa( previous,childClass )
                self.children{ind} = [];
            else
                self.children{end+1} = child;
            end
            
            % now update the child
            for j = 1:numel( child )
                parentClass = class( self );
                [previous,ind] = child(j).getParent( parentClass );
                
                if ~isempty( previous ) && isa( previous,parentClass )
                    child(j).parent{ind}(end+1) = self;
                elseif isempty( previous ) && isa( previous,parentClass )
                    child(j).parent{ind} = [];
                else
                    child(j).parent{end+1} = self;
                end
            end
        end
        
        
        function removeChild( self,type,varargin )
            % removeChild( self,type,(indices) )
            %
            % remove a child of class "type" and optionally specific
            % indices of that child class

            % get children of the given type
            [child,idx] = self.getChild( type );
            if ~isempty( idx ) && max( idx ) >= 1
                if nargin == 3
                    indices = varargin{1};
                else
                    indices = 1:numel( child );
                end
            
                % first remove the parent from the child
                for j = indices
                    [~,parentID] = child(j).getParent( class(self) );
                    if ~isempty( parentID ) % if it is, "self" is not a parent
                        child(j).parent{parentID} = [];
                    end
                end
                
                % now remove the child from the parent
                if nargin == 3 || numel( indices ) == numel( child )
                    self.children{idx}(indices) = [];
                else
                    allInd = true( 1,numel(self.children ) );
                    allInd(idx) = false;
                    self.children = self.children(allInd); % remove all instances of that child
                end
                
                % clean up the parent & child
                self.children = self.children(~cellfun( @isempty,self.children ));
                for j = 1:numel( child )
                    child(j).parent = child(j).parent(~cellfun( @isempty,child(j).parent ));
                end
            end 
        end
       
        
        function [child,indices] = getChild( self,type,varargin)
            % [child,indices] = getChild( self,type,(index) )
            %
            % pull out the children of a given type, and optionally,
            % only those with a specific index
            child = [];
            indices = find( cellfun( @(c)(isa( c,type )),self.children ) );
            if ~isempty( indices )
                child = self.children{indices};
                if nargin == 3
                    child = child(varargin{1});
                end
            end 
        end   
           
        
        function addParent( self,parent )
            % addParent( self,parent )
            %
            % update the parents of an object, and update that parent to
            % reference the new child

            % find the previous parents of the given type
            previous = cell2mat( cellfun( @(c)(isa( c,class( parent ) )),...
                self.parent,'un',0 ) );
            if  isempty( previous ) || max( previous ) == 0
                self.parent{end+1} = parent;
            else
                self.parent{previous}(end+1) = parent;
            end

            % now update the parent 
            previous = cell2mat( cellfun( @(c)(isa( c,class( self ) )),...
                parent.children,'un',0 ) );
            if isempty( previous ) || max( previous ) == 0
                parent.children{end+1} = self;
            else
                parent.children{previous}(end+1) = self;
            end
        end
        
        
        function removeParent( self,type,varargin )
            % removeParent( self,type,(indices) )
            %
            % remove a parent of class "type", and optionally, the specific
            % indices of that parent class

            % find the parents of class "type"
            [padre,idx] = self.getParent( type );
            if ~isempty( idx ) && max( idx ) == 1
                if nargin == 3
                    indices = varargin{1};
                else
                    indices = 1:numel( padre );
                end
            
                % first remove the child from the parent
                for j = indices
                    [child,childID] = padre(j).getChild( class( self ) );
                    thisChild = find( ismember( child,self ) );                    
                    padre(j).children{childID}(thisChild) = [];
                end
                clear child 
                
                % now remove the parent from the child
                if nargin == 3 || numel( indices ) == numel( padre )
                    self.parent{idx}(indices) = [];
                else
                    allInd = true( 1,numel( self.parent ) );
                    allInd(idx) = false;
                    self.parent = self.parent(allInd); % remove all instances of this parent
                end
                
                % clean up the parent & child
                for j = 1:numel( padre )
                    padre(j).children = padre(j).children(~cellfun( @isempty,padre(j).children ));
                end
                self.parent = self.parent(~cellfun( @isempty,self.parent ));
            end 
        end
        
        
        function [parent,indices] = getParent( self,type,varargin )
            % [parent,indices] = getParent( self,type,(index) )
            %
            % pull out the parents of a given type, and optionally,
            % those with a specific index
            parent = [];
            indices = find( cellfun( @(c)(isa( c,type )),self.parent ) );
            if ~isempty( indices )
                parent = self.parent{indices};
                if nargin == 3
                    parent = parent(varargin{1});
                end
            end
        end
        
        
        function sibling = getSibling( self,type,parentType )
            % sibilng = getSibling( self,type,parentType )
            %
            % pull out the siblings of a given type (i.e. those that share
            % a common parent)
            %
            % loop over the parents contained within the current object.
            % For each, find all children (excluding the current object) of
            % the class specified by "type". Only get children
            % of the parent of type "parentType" as a second input.
                        
            sibling = [];
            parents = self.getParent( parentType );
            if isempty( parents )
                error( ['No parents of type ',parentType,' found'] );
            end

            % loop over parents
            for j = 1:numel( parents )
                child = parents(j).getChild( type ); % gets the siblings
                if ~isempty( child )
                    child(child == self) = []; % removes redundancy
                    sibling = [sibling, child];
                end
            end
        end


        function partner = getPartner( self,type,childType )
            % partner = getPartner( self,type,childType )
            %
            % pull out the partners of a given type (i.e. those that share 
            % a common child)
            %
            % loop over the children contained within the current object. 
            % For each child of type "childType", find all parents
            % (excluding the current object), and extract those of type "type"
            
            partner = [];
            child = self.getChild( childType );
            if isempty( child )
                error( ['No children of type ',childType,' found'] )
            end

            % loop over the children
            for j = 1:numel( child )
                padre = child(j).getParent( type );
                if ~isempty( padre )
                    padre(padre == self) = []; % removes self reference
                    partner = [partner, padre];
                end
            end
        end

    
        function deleteSelf( self )
            % deleteSelf( self )
            %
            % delete this object from memory. Remove from parents and 
            % remove all children of this object
            for j = 1:numel( self.parent )
                self.removeParent( class( self.parent{1} ) );
            end
            for j = 1:numel( self.children )
                self.removeChild( class( self.children{1} ) );
            end
            delete( self );
        end
        
    end % methods
    
end