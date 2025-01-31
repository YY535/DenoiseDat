function DenoiseDatYY(datFile,session,varargin)

% DenoiseDat - remove the first PCA ('noise') component from your .dat file
% 
% This function computes PCA (each shank is treated separately) over a
% short (preferably noisy) sample of the data (from the period "baseline", 
% user-selected to include noisy data). Subsequently, the part of the 
% signal that is explained by this "noisy" component is removed, resulting
% in a denoised ('cleaned') .dat file, which should improve spike sorting.
% Note: the selection of a "baseline" period containing the problematic noise
% (physical or electrical noise) is crucial for the denoising to work. 
%
% Note that this manipulation will also remove any real signal shared among
% channels, such as important oscillations. It is therefore recommended
% to NOT apply this de-noising to the original dat files as this will 
% permanently modify them and delete this signal. Instead, it is recommended
% that you preserve your original files in their subsession folders, and 
% you apply this denoising step to the .dat file only after the LFP file
% has been created (the denoised .dat file will be relatively flat, missing
% a lot of oscillation information of interest for the .lfp). 
%
%  USAGE
%
%    DenoiseDat(session, datFile, <options>)
%
%  INPUTS
%  
%    datFile        full path towards the file to be denoised. 
%                   Make sure original files are backed up!
%    session        Cell Explorer format variable containing parameters about 
%                   the recording (sampling rate, shank configuration, number 
%                   of channels)
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties         Values
%    -------------------------------------------------------------------------
%   'SSD_path'           path to SSD. If provided, .dat file will be copied
%                        and edited there for faster performance, and 
%                        afterwards copied back, overwriting the original file.
%   'baseline'           interval (in seconds) that will be used to compute
%                        the PCA. It is recommended that this includes a 
%                        noisy period (default = whole recording)
%   'sampleDuration'     duration (in seconds) of the random sample that 
%                        will be selected from the "baseline" period to 
%                        compute the PCA (default = 20). Note that longer 
%                        periods increase memory demands. 
%   'secondsPerChunk'    duration (in seconds) of a chunk of data (default = 10). 
%                        After computing the noise component, noise is 
%                        progressively removed from the .dat file in chunks. 
%                        Note that larger chunks increase memory demands. 
%   'FlatThreshold'     thereshold to pickup flat components: defualt =
%                       10 (std of the spatial loading variance)
%   'rejectChannels'     channels (base1) that should be excluded from denoising
%   'verbose'            Display progress messages (default = true)
%   'groupshanks'        groups of shanks to denoise jointly. pls don't
%                        provide overlapping in groups. {[shk1,shk2], [shk3 shk4 shk5]}
%                       it usually works better if you group shanks in the
%                       same region, recommend to have channel number >30. 
%   'silentperiod'      the periods denoising skips, e.g. sleeping periods.
%                       default [], then denoise all the periods. 
%                       shape: nperiods x 2, [starts, ends] in s.
%   'fitforwhole'       one component for the whole data. default: false
%   'isoverwrite'       if overwrite the dat file. default: false
%    =========================================================================
%
%  EXAMPLE 
%  
%    LFPfromDat; % Produce LFP before denoising the .dat file to preserve slow oscillations affecting all channels (e.g. delta waves)
%    DenoiseDat(session,datFile,'SSD_path','D:\'); % remove first PCA 'noise' component from .dat file to improve spike sorting 
%    kilosortFolder = KiloSortWrapper('SSD_path','D:\','rejectchannels',excludeChannels,'datFilename',SSD_file); % spike sort
%
% IO: Ralitsa Todorova and Weiwei Chen 2024, Method: Weiwei 2024
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%-------------------------------------------------------------------------

p = inputParser;
addParameter(p, 'SSD_path', [], @isfolder);
addParameter(p, 'baseline', [], @isdmatrix); % interval to use as basis for PCA
addParameter(p, 'sampleDuration', 20, @isdscalar); % cap duration of the noise random sample to be used to compute PCA
addParameter(p, 'secondsPerChunk', 120, @isdscalar);% every 10s. well, it's ok
addParameter(p, 'FlatThreshold', 10, @isdscalar);
addParameter(p, 'rejectChannels', [], @isivector);
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'groupshanks', {}, @iscell);
addParameter(p, 'silentperiod', [], @ismatrix);
addParameter(p, 'fitforwhole', false, @islogical);
addParameter(p, 'isoverwrite', false, @islogical);

parse(p, varargin{:});

SSD_path = p.Results.SSD_path;
baseline = p.Results.baseline;
sampleDuration = p.Results.sampleDuration;
secondsPerChunk = p.Results.secondsPerChunk;
flatthr = p.Results.FlatThreshold; 
rejectChannels = p.Results.rejectChannels; %!
verbose = p.Results.verbose;
groupshanks = p.Results.groupshanks;
silentperiod = p.Results.silentperiod * session.extracellular.sr;
fitforwhole = p.Results.fitforwhole;
isoverwrite = p.Results.isoverwrite;
%% Load information from the session variable
nChannels = session.extracellular.nChannels;
samplingRate = session.extracellular.sr;
shanks = session.extracellular.spikeGroups.channels; nShanks = length(shanks);
for k=1:nShanks, shanks{k}(ismember(shanks{k},rejectChannels)) = []; end % remove any channels from this process
% recollect the channels and bin the shanks for jointly denoising
if ~isempty(groupshanks)
    orig_shanks = shanks;
    ngroup = length(groupshanks);
    shanks = cell(1,ngroup);
    shanks_ch = cell(1,ngroup);
    channel_shank = cell(1,nShanks);
    for k = 1:nShanks
        channel_shank{k} = k*ones(size(orig_shanks{k}));
    end

    for k=1:ngroup
        shanks{k} = cell2mat(orig_shanks(groupshanks{k}));

        shanks_ch{k} = cell2mat(channel_shank(groupshanks{k}));
    end
    nShanks = length(shanks);
end
basepath = session.general.basePath;
basename = basenameFromBasepath(basepath);

% Create .EMG.dat file for removed noise activity and the .EMGcomp.mat for
% nosie spactial loading. 
EMG_noisefile = fullfile(SSD_path,[basename '.EMG.dat']);
EMG_comps = fullfile(SSD_path,[basename '.EMGcomp.mat']);

% the original .dat file could be restored by adding A*s back. 

% m = memmapfile(datFile, 'Format','int16','Writable',true);
if ~exist(datFile,'file'), error([datFile ' file does not exist.']); end

if isempty(SSD_path), SSD = false; else, SSD = true; end
if isoverwrite
    SSD_file = datFile;
    fprintf('\n Overwriting the original dat file: \n %s\n', SSD_file)
else
    SSD_file = fullfile(SSD_path,[basename 'D.dat']);copyfile(datFile,SSD_file);
end
m = memmapfile(SSD_file, 'Format','int16','Writable',true);

data = reshape(m.data,nChannels,[]);
nSamples = size(data,2);



% Compute the PCA shank by shank
% notice, for shanks<32 usually recommend to compute use multiple shanks. 
V = cell(nShanks,1); vs = cell(nShanks,1);
% EMGcomp = nan(nShanks,1);
opf_norm = @(x)x/sqrt((x'*x));
% 
if fitforwhole
% fit a simple model to the whole period. 
    if verbose, disp([datestr(clock) ': Computing PCA over noise periods in ' basepath '...']); end

    % Select timestamp indices that will be taken as baseline to compute PCA:
    if isempty(baseline), baseline = [0 nSamples/samplingRate]; end % use the whole recording as a "baseline"
    baselineSamples = round(baseline*samplingRate);
    baselineSamples(diff(baselineSamples,[],2)<=0,:) = [];
    if sum(diff(baselineSamples,[],2))>samplingRate*sampleDuration % baseline is larger than the required number of indices
        % select a random sample of indices within the baseline period
        idx = sort(randperm(round(sum(diff(baselineSamples,[],2))),round(samplingRate*sampleDuration)))'; % select random indices among all possible indices
        idx = round(Unshift(idx,baselineSamples)); % shift them to the corresponding indices within the baseline
    else
        if isdvector(baselineSamples)
            idx = baselineSamples(1):baselineSamples(2);
        else
            idx = linspaceVector(baselineSamples(:,1),baselineSamples(:,2)); % all possible indices within the baseline period
        end
    end
    if ~isempty(silentperiod)
        for k = 1:size(silentperiod,1)
            idx(idx>=silentperiod(k,1)&idx<silentperiod(k,2))=[];
        end
    end
    if isempty(idx)
        warning('Skipping all the periods.')
        return
    end
    baselineData = data(:,idx);
    % find flat components with PCA 
    Vprojection = zeros(nChannels,nShanks);
    for k=1:nShanks
        channels = shanks{k};
        %%
        u = nan(length(channels),5);
        for kperm = 1:5
            shankData = baselineData(channels,ceil(length(idx)*rand(1e5,1)));
            % new version: add flatness control 
            [tmp_u,~,~] = svd(cov(double(shankData')));
            flat_u = flatness(tmp_u);
            [v_k, tmp_k] = max(flat_u);
            % tmp_k = find(flat_u>(flatthr*std(tmp_u)));
            if v_k>=(flatthr*std(tmp_u))
                u(:,kperm) = tmp_u(:,tmp_k);
            end
        end
        V{k} = u;
        vs{k} =opf_norm(nanmean(bsxfun(@times,u,sign(median(u))),2));
        Vprojection(channels,k) = vs{k};
    end
    % out.EMGcomp = EMGcomp;
    out.V = V;
    out.vs = vs;
    out.shanks = shanks;
    out.baselinesamples = single(idx);
    save(EMG_comps,'out')
    clear out baselineData
end

if verbose, disp([datestr(clock) ': Removing noise components in ' datFile ' file...']); end
%%
% Modify the .dat file in chunks:
chunkSize = ceil(samplingRate*secondsPerChunk);
indicesChunk1 = (1:chunkSize*nChannels)'; % these are the .dat file indices for the first chunk
nChunks = ceil(nSamples/chunkSize);

beginningend = bsxfun(@plus, [1:nChunks]'*chunkSize, [1-chunkSize 0]);
beginningend(end) = min(beginningend(end),nSamples);

% skip large sleeping periods. 
useprd = [1:nChunks]';
slpperiods = cell(nChunks,1);
if ~isempty(silentperiod)
    for k=1:size(silentperiod,1)
        useprd(beginningend(:,1)>=silentperiod(k,1)&beginningend(:,2)<=silentperiod(k,2))=0;
        bgtmp = find(beginningend(:,1)<=silentperiod(k,1)&beginningend(:,2)>=silentperiod(k,1));
        slpperiods{bgtmp} = [slpperiods{bgtmp}; [silentperiod(k,1) min(silentperiod(k,2), beginningend(bgtmp,2))]];
        endtmp = find(beginningend(:,1)<=silentperiod(k,2)&beginningend(:,2)>=silentperiod(k,2));
        slpperiods{endtmp} = [slpperiods{endtmp}; [max(silentperiod(k,1), beginningend(endtmp,1)), silentperiod(k,2)]];

    end
end
EMG_s = zeros(nShanks,nSamples,'int16');


for chunk = 1:nChunks
    if ~useprd(chunk) % skip silent periods
        continue
    end
    tmp_ = beginningend(chunk,1):beginningend(chunk,2);
    slpperiods{chunk} = BinPeriod(slpperiods{chunk},1);
    dataChunk = (data(:,tmp_));

    indicesChunk = indicesChunk1+(chunkSize*nChannels)*(chunk-1); if indicesChunk(end)>nSamples*nChannels, indicesChunk(indicesChunk>nSamples*nChannels) = []; end
    noise = zeros(size(dataChunk),'int16');
    for k=1:nShanks
        tmp_ch = unique(shanks_ch{k});
        mdataChunk = double(dataChunk(shanks{k},:));% 3ch median filter+ highpassed, to get the surrogate noisy data. 
        for kk=1:length(tmp_ch)
            mdataChunk(shanks_ch{k}==tmp_ch(kk),:) = medfilt1(mdataChunk(shanks_ch{k}==tmp_ch(kk),:),3);
        end
        mdataChunk = ButFilter(mdataChunk',4,80/(2e4/2),'high');% phase differenrce removed
        if fitforwhole
            tmp_u = vs{k};
        else
            [u,~,~] = svd(cov(mdataChunk));% compute the EMG noise from the surrogate data
            flat_u = flatness(u);
            [~,use_u ] = max(flat_u);
            vs{k,chunk} = u(:,use_u)*sign(sum(u(:,use_u)));
            tmp_u = vs{k,chunk} ;
        end
        projection =  reshape(detrend(mdataChunk*tmp_u),1,[]);
        for kp = 1:size(slpperiods{chunk},1)
            tmp_slp = slpperiods{chunk}(kp,:)-beginningend(chunk,1)+1;
            dtmp_slp = diff(tmp_slp);
            [~,idxbg] = min(abs(projection(tmp_slp(1)+(0:min(200, dtmp_slp)))));
            [~,idxend] = min(abs(projection(tmp_slp(2)+(max(-199, -dtmp_slp):0))));
            projection((tmp_slp(1)-1+idxbg):(tmp_slp(2)-200+idxend))=0;
        end
        % project back to channels:
        noise(shanks{k},:) = int16(vs{k}*projection);
        EMG_s(k,tmp_) = projection;
    end
    m.Data(indicesChunk) =   m.Data(indicesChunk)-noise(:);
end
if ~fitforwhole
out.V = [];
out.vs = vs;
out.shanks = shanks;
out.shanks_ch = shanks_ch;
out.beginningend = beginningend;
out.slpperiods = slpperiods;
save(EMG_comps,'out')
clear out 
end

fileID = fopen(EMG_noisefile,'w');
fwrite(fileID, EMG_s(:),'int16');
fclose(fileID);
disp([datestr(clock) ': Noise removed!']);

clear m

