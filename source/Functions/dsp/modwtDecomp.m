function [wave,power,frequencies] = modwtDecomp( data,Fs,wtype,wlevel )
% [wave,power,frequencies] = modwtDecomp( data,Fs,wtype,wlevel )
%
% Computes the maximal-overlap discrete wavelet transform (modwt) over each
% column in the n x m data matrix "data" using the wavelet type specified
% by "wtype", and the wavelet level specified by "wlevel".
%
%                   >>> INPUTS >>>
% data: 
%   n x m data matrix in column-major format
% Fs:
%   the sampling rate
% wtype:
%   string specifying the wavelet type 
% wlevel:
%   the final wavelet level to decompose 
%
%                   <<< OUTPUTS <<<
% wave:
%   n x wlevel x m tensor of wavelet coefficients from the modwt
% power:
%   |wave|^2
% frequencies:
%   wlevel+1 x 2 matrix specifying the frequency bandwidth in each level
% 
% by JMS, 5/20/2016

% maximal-overlap discrete wavelet transform (MODWT)
[n,m] = size( data );
wave = zeros( wlevel+1,n,m );
disp('Computing wavelet');
for j = 1:m
    fprintf(' . ');
    wave(:,:,j) = modwt( data(:,j),wtype,wlevel );
end
power = abs( wave ).^2; % compute wavelet power spectrum
fprintf('\nWavelet calculated\n');

% compute the frequency bands represented by the transform
frequencies = zeros( wlevel+1,2 ); 
levels = 1:wlevel+1;
frequencies(:,1) = 0.5 * Fs ./ (2.^levels); % 1/2 (2^-j * Fs)
frequencies(:,2) = 0.5 * Fs ./ (2.^(levels-1));
frequencies(end,1) = 0;
frequencies(end,2) = Fs/2/(2^wlevel);

end