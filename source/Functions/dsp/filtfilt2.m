% --------------  [fd] = filtfilt2(d,hp,lp,Fs,varargin) ----------------------
%
% Like "filtfilt", this function performs a forward and reverse filter to remove
% the phase lag induced by a single forward filter.  Unlike filtfilt, which
% uses [b,a] methodology, this function utilizes [z,p,k] coefficients which
% are more stable at extreme low/high filter cutoffs.
%
%               >>> INPUTS >>>
%   REQUIRED:
%       d = data to be filtered
%       hp = highpass cutoff (0 = no highpass)
%       lp = lowpass cutoff (0 = no lowpass)
%       Fs = sampling rate
%   
%   OPTIONAL:
%       hpR = highpass roll off (default 6)
%       lpR = lowpass roll off (default 6)
%       type = 'butter','cheby1','ellip' (default 'butter')
%       plotting = 1 or 0
%           plots data and filtered data; plots filter design
%           (default 0)
%       mat = 'row' or 'col' 
%           specifies if data is in row or column vector format
%           (default 'col')
%
%               <<< OUTPUTS <<<
%   fd = filtered data
%
%
%   Example:
%       Fs = 1000;
%       si = 1/Fs;
%       time = 0:si:1-si;
%       noise = cos((2*pi)*500*time);
%       sig = 2*cos((2*pi)*10*time)+noise;
%       hp = 1;
%       lp = 100;
%       filtSig = filtfilt2(sig,hp,lp,Fs);
%       plot(time,sig,time,filtSig,'r','linewidth',2);
%
% By JMS, 2/12/2015
%------------------------------------------------------------

function [fd] = filtfilt2(d,hp,lp,Fs,varargin)

% extract varargins and defaults
if nargin<9;mat='col';
else mat=varargin{5};end
if nargin<8;plotting=0;
else plotting=varargin{4};end
if nargin<7;type='butter';
else type=varargin{3};end
if nargin<6;lpR=6;
else lpR=varargin{2};end
if nargin<5;hpR=6;
else hpR=varargin{1};end

if plotting>0;plotting=1;end
if ~strcmpi(type,{'butter','cheby1','ellip'});type='butter';end
if lpR<0;lpR=6;end
if hpR<0;hpR=6;end
if isempty(plotting);plotting=1;end
if isempty(type);type='butter';end
if isempty(hpR);hpR=6;end
if isempty(lpR);lpR=6;end

% design filter 
if strcmpi(type,'butter')
    if lp && hp
       [z,p,k] = butter(hpR,[hp lp]/(Fs/2));
    elseif lp
        [z,p,k] = butter(lpR,lp/(Fs/2));
    elseif hp
        [z,p,k] = butter(hpR,hp/(Fs/2),'high');
    end
elseif strcmpi(type,'cheby1')
    if lp && hp
       [z,p,k] = cheby1(hpR,[hp lp]/(Fs/2));
    elseif lp
        [z,p,k] = cheby1(lpR,lp/(Fs/2));
    elseif hp
        [z,p,k] = cheby1(hpR,hp/(Fs/2),'high');
    end
elseif strcmpi(type,'ellip')
     if lp && hp
       [z,p,k] = ellip(hpR,[hp lp]/(Fs/2));
    elseif lp
        [z,p,k] = ellip(lpR,lp/(Fs/2));
    elseif hp
        [z,p,k] = ellip(hpR,hp/(Fs/2),'high');
     end
end

% transpose data if entered as row-vector format
if min(size(d)) > 1
    if strcmpi(mat,'row');
        d = d'; % transpose
    end
elseif size(d,1) == 1
    d = d'; % make column vector
end

% filter data and reverse filter to remove phase shift
[s,g] = zp2sos(z,p,k);
Hd = dfilt.df1sos(s,g);
fd = filter(Hd,d);
fd = flipud(fd);
fd = filter(Hd,fd);
fd = flipud(fd);

% plot data, filtered data, and filter design
if plotting
    figure; hold on;
    plot(d,'k');
    plot(fd,'r');
    legend({'raw','filtered'});
    fvtool(s,g,'Analysis','freq');
end

end
    