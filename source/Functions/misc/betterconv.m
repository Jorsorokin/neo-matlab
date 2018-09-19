function Y = betterconv( X,kernel )
    % Y = betterconv( X,kernel )
    %
    % pads the data in "X" with zeros the length of the kernel before
    % convolving to avoid edge artifacts. also ensures the kernel is
    % sum-normalized. X will be convolved along the first dimension (rows)

    kernel = kernel/sum( kernel );
    padSize = length( kernel );
    [n,m,p] = size( X );
    pad = zeros( padSize,m,p,class( X ) );

    Y = convn( [pad;X;pad],kernel,'same' );
    Y = Y(padSize+1:n+padSize,:,:);
end