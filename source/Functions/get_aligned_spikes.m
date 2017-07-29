eegch = 'I1';
spch = [1:4];
win = 4; % starting window time
w = 2; % pre/post time window
runstats = 1; % for plotting and doing wilcoxon signed-ranks tests
align_manually = 0 % to align manually...if not, uses "leadtime" from the R structs
separate_shamlaser = 1; % to plot statistics separately for sham/laser

T = dataset;
files = uipickfiles;
% LOOP OVER FILES
for f = 1:numel(files)
    try
    
        clear sPre sPost bPre bPost sptm bt bs sRate bRate ISI IBI R S mu mufs spsnip snip spikes sp ref good bad
    
        % load the eeg results, get spikes
        cd(files{f});
        R = struct2array(load([files{f},'/eegResults.mat'])); % eeg
        if isfield(R.trialInfo,'chanlocs')
            chan = find(strcmp(R.trialInfo.chanlocs,eegch));
        else
            chan = 1;
        end
        close;

        % get the eeg data
        if isfield(R,'validSeizures')
            if isfield(R.validSeizures,'laserConcat') || isfield(R.validSeizures,'shamConcat')
                vS = R.validSeizues.shamIndex;
                vL = R.validSeizures.laserIndex;
            else               
                vS = boolean(R.validSeizures.shamValid);
                vL = boolean(R.validSeizures.laserValid);
            end
        else
            vS = boolean(ones(sum(R.allevents(:,end)==0),1)); 
            vL = boolean(ones(sum(R.allevents(:,end)==1),1));
        end

        % sampling rate, and # points for window
        if isfield(R,'trailInfo')
            Fs = R.trailInfo.Fs;
        else
            Fs = R.trialInfo.Fs;
        end

        % extract eeg segments, starting/ending times of the segments
        if min(R.allevents(:,end))==1
            if isfield(R.validSeizures,'laserConcat')
                eeg = cellfun(@(c)(c(:,chan)),R.validSeizures.laserConcat,'un',0);
                n = numel(eeg{1});
                eeg = reshape(cell2mat(eeg),n,numel(eeg));
            else
                eeg = R.laserEvents{chan}(vL,:)';
                ind = horzcat(R.times.startIndex{1}(vL), R.times.endIndex{1}(vL));
            end
            vS = [];
        elseif max(R.allevents(:,end))==0
            if isfield(R.validSeizures,'shamConcat')
                eeg = cellfun(@(c)(c(:,chan)),R.validSeizures.shamConcat,'un',0);
                n = numel(eeg{1});
                eeg = reshape(cell2mat(eeg),n,numel(eeg));
            else
                eeg = R.shamEvents{chan}(vS,:)';
                ind = horzcat(R.times.startIndex{2}(vS), R.times.endIndex{2}(vS));
            end
            vL = [];
        else    
            if isfield(R.validSeizures,'laserConcat')
                L = cellfun(@(c)(c(:,chan)),R.validSeizures.laserConcat,'un',0);
                S = cellfun(@(c)(c(:,chan)),R.validSeizures.shamConcat,'un',0);
                n = numel(eeg{1});
                
                L = reshape(cell2mat(L),n,numel(L));
                S = reshape(cell2mat(S),n,numel(S));
                eeg = [L,S];
            else
                eeg = [R.laserEvents{chan}(vL,:);R.shamEvents{chan}(vS,:)]';
                ind = horzcat([R.times.startIndex{1}(vL); R.times.startIndex{2}(vS)], [R.times.endIndex{1}(vL); R.times.endIndex{2}(vS)]);
            end
        end
        ind = ind/Fs; % convert to seconds
        key = logical(vertcat(ones(sum(vL),1), zeros(sum(vS),1)));

        % pull out multi-unit data for the given start/stop times
        [mu,mufs] = open_tdt_data(R.fullfile,spch,[],[],'MUA1');
        spikes = zeros(numel(ind(1,1)*mufs:ind(1,2)*mufs),size(eeg,2),size(mu,2));
        for c = 1:size(mu,2)
            for j = 1:size(ind,1)
                x = mu(ind(j,1)*mufs:ind(j,2)*mufs,c);
                if numel(x) > size(spikes,1)
                    x = x(1:end-1);
                elseif numel(x) < size(spikes,1)
                    x = [x;0];
                end   
                spikes(:,j,c) = x;
                clear x
            end
        end
        
        % choose the best spike channel, and reference against the worst to
        % eliminate common-mode noise
        subplot(2,1,1);
        multisignalplot(squeeze(spikes(:,1,:)),mufs);
        subplot(2,1,2);
        multisignalplot(squeeze(spikes(:,size(eeg,2),:)),mufs);
        good = input('Use which ch? ');
        bad = input('Reference which ch? ');
        sp = squeeze(spikes(:,:,good));
        if bad ~= 0
            ref = squeeze(spikes(:,:,bad));
            sp = sp - ref;
        end
        close
        

        % get the aligned seizures and spikes
        %========================================================
        if align_manually == 1
            if isfield(R.times,'alignedtimes')
                [~,~,snip,spsnip] = alignszr(win,eeg,Fs,sp,mufs,R.times.alignedtimes(:,1),R.times.alignedtimes(:,2));
            else
                [start,stop,snip,spsnip] = alignszr(win,eeg,Fs,sp,mufs);
                R.times.alignedtimes = [start,stop];
                results = R;
                save('eegResults.mat','-append','results');
                clear results
            end
        else
            % else use the "leadtime" from the "R" structure to pull out
            % eeg and spike snips
            start = R.times.leadTime/Fs - win;
            stop = R.times.leadTime/Fs + win;
            snip = eeg(start*Fs:stop*Fs,:); 
            spsnip = sp(start*mufs:stop*mufs,:);
            figure; multisignalplot(snip,Fs);
            figure; multisignalplot(spsnip,mufs);
        end
            
        %========================================================
        
        % determine which ones to eliminate
        bad = input('Exclude which segments? ');
        if bad ~= 0
            snip(:,bad) = [];
            spsnip(:,bad) = [];
            key(bad) = [];
        end
        close all

        % Re-find spikes and bursts on the aligned spike snips
        sptm = SpikeDetector(spsnip,mufs,'thresh',2.2,'artifact',0.4,'plottrace',1);
        sptm = sptm./mufs;

        % call "findbursts" to calculate burst times
        %=========================
        [bs,bt] = findbursts(sptm*mufs,0.01,5,1/mufs,0.06,3,0.01);
        %=========================

        % plot rasters of spikes and bursts
        figure;
        subplot(2,1,1);
        PlotRasters(sptm,win,w,w,[],0); % spikes
        subplot(2,1,2);
        PlotRasters(bt,win,w,w,[],0); % bursts
        pause(4);
        close all; 
        
        % spike/burst rates for pre/post
        sPre = (sptm >= win-w & sptm < win);
        sPost = (sptm >= win & sptm <= win+w);
        bPre = (bt >= win-w & bt < win);
        bPost = (bt >= win & bt <= win+w);

        sRate(:,1) = sum(sPre)' ./ w;
        sRate(:,2) = sum(sPost)' ./ w;
        bRate(:,1) = sum(bPre)' ./ w;
        bRate(:,2) = sum(bPost)' ./ w;

        % ISI/IBI
        for j = 1:size(bt,2)
            x = diff(sptm(sPre(:,j),j));
            y = diff(sptm(sPost(:,j),j));
            x(x<=0) = [];
            y(y<=0) = [];
            ISI(j,1) = nanmean(x);
            ISI(j,2) = nanmean(y);

            x = diff(bt(bPre(:,j),j));
            y = diff(bt(bPost(:,j),j));
            x(x<=0) = [];
            y(y<=0) = [];
            IBI(j,1) = nanmean(x);
            IBI(j,2) = nanmean(y);
        end

        % store into dataset T
        T.Folder{f,1} = files{f};
        T.Virus{f,1} = R.trialInfo.subStrain;
        T.Cage{f,1} = R.trialInfo.cage;
        T.Date{f,1} = R.trialInfo.date;
        T.Trial{f,1} = R.trialInfo.trialNum;
        T.Type{f,1} = R.trialInfo.trigger;
        T.StimKey{f,1} = key; 
        T.eegChan{f,1} = eegch; 
        T.spChan(f,1) = good;
        T.Removed{f,1} = bad;
        T.EEG{f,1} = snip;
        T.Spikes{f,1} = spsnip;
        T.spTimes{f,1} = sptm;
        T.btTimes{f,1} = bt;
        T.btStruct{f,1} = bs;
        T.spRate{f,1} = sRate;
        T.ISI{f,1} = ISI;
        T.btRate{f,1} = bRate;
        T.IBI{f,1} = IBI;
        T.Index(f,1) = f;
    catch
        disp(lasterr);
    end
end

%% Plot some stuff
clear spRate btRate ISI IBI spVar btVar
i = 1;
for j = 1:size(T,1)
    if T.spChan(j) > 0
        i = i+1;
        LspRate(i,:) = nanmean(T.spRate{j}(T.StimKey{j},:));
        LbtRate(i,:) = nanmean(T.btRate{j}(T.StimKey{j},:));
        LspVar(i,:) = nanstd(T.spRate{j}(T.StimKey{j},:));
        LbtVar(i,:) = nanstd(T.btRate{j}(T.StimKey{j},:));
        LISI(i,:) = nanmean(T.ISI{j}(T.StimKey{j},:));
        LIBI(i,:) = nanmean(T.IBI{j}(T.StimKey{j},:));

        SspRate(i,:) = nanmean(T.spRate{j}(~T.StimKey{j},:));
        SbtRate(i,:) = nanmean(T.btRate{j}(~T.StimKey{j},:));
        SspVar(i,:) = nanstd(T.spRate{j}(~T.StimKey{j},:));
        SbtVar(i,:) = nanstd(T.btRate{j}(~T.StimKey{j},:));
        SISI(i,:) = nanmean(T.ISI{j}(~T.StimKey{j},:));
        SIBI(i,:) = nanmean(T.IBI{j}(~T.StimKey{j},:)); 
    end
end
       
if separate_stimlaser == 1
    pSp(1,1) = signrank(LspRate(:,1),LspRate(:,2));
    pSp(2,1) = signrank(LspVar(:,1),LspVar(:,2));
    pSp(3,1) = signrank(LISI(:,1),LISI(:,2));
    pSp(1,2) = signrank(SspRate(:,1),SspRate(:,2));
    pSp(2,2) = signrank(SspVar(:,1),SspVar(:,2));
    pSp(3,2) = signrank(SISI(:,1),SISI(:,2));
    
    pBt(1,1) = signrank(LbtRate(:,1),LbtRate(:,2));
    pBt(2,1) = signrank(LbtVar(:,1),LbtVar(:,2));
    pBt(3,1) = signrank(LIBI(:,1),LIBI(:,2));
    pBt(1,2) = signrank(SbtRate(:,1),SbtRate(:,2));
    pBt(2,2) = signrank(SbtVar(:,1),SbtVar(:,2));
    pBt(3,2) = signrank(SIBI(:,1),SIBI(:,2));
end
    
