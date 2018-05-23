function features = compute_wavelet_features( X,level,type )
    features = [];
    for j = 1:size( X,3 )
        dec = mdwtdec( 'c',X(:,:,j),level,type );
        coeffs = cellfun( @(x)(x'),dec.cd,'un',0 );
        features(:,:,j) = [coeffs{:}];
    end
end
            
        