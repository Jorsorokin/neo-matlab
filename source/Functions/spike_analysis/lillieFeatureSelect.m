function [scores,hits] = lillieFeatureSelect( features )

% pre allocate
scores = zeros( size( features,2 ),1 );
hits = boolean( scores );

% loop over the feature set
for j = 1:size( features,2 )
    [~,~,scores(j)] = lillietest( features(:,j) );
end

% get the noise level
noise = std( scores );

% find those with large scores
hits( scores >= 2*noise ) = true;

end