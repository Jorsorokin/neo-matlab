function [split,V] = check_voltage_distribution( X )
    % [split,V] = check_voltage_distribution( X )
    %
    % checks whether the join distribution of the minimum voltage on
    % each channel across all spikes is unimodal (i.e. should follow a
    % chi-squared distribution). If any relevant channel shows
    % non-unimodality, then split becomes true.
    %
    % Inputs:
    %   X - an n x m x c matrix, with n = points, m = spikes (observations)
    %       and c = channels
    %
    % Outputs:
    %   split - a boolean of whether or not the voltage is bimodal or not.
    %           If true, it means the voltage distribution is not unimodal,
    %           and so likely this set of spike waveforms should be split
    %
    %   V - the minimum voltages on relevant channels
    %
    % By Jordan Sorokin, 7/16/18

    v = var( squeeze( mean( X,2 ) ) ); 
    V = squeeze( min( X(:,:,v > mad(v)) ) );
    
    % test for unimodal, multi-variate gaussian via chi2 test
    [~,pval,cramerV] = chi2normal( V );
    split = pval < 0.001 & cramerV > 0.3;
end