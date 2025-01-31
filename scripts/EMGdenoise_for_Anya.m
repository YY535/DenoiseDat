%% example: cleaning multiple session from one animal without sleeping period
% for now pls addpath to ayadata4: 
% Weiwei\drafts\util 
% and 
% Weiwei\drafts\AnalysisValidations

database = 'X:\AGRP\EphysAgRP\MCh7\';
sessionnames = {'day8','day10','day11'};

savebase = 'D:\WorkingDir\MCh7\';% replace with your working dir

isoverwrite_dat = true;% if you want to overwrite the orignial dat file. defualt: false
% otherwise the file would be saved in SSD_path
Sleepprd = [];
for k = 1:length(sessionnames)
    %%
    sessionname = sessionnames{k};
    SSD_path = [savebase,sessionname];
    if ~exist(SSD_path,"dir")
        mkdir(SSD_path)
    end
    filebase = [database,'\',sessionname,'\'];
    filename = [filebase,'\',sessionname];
    
    load([filename, '.session.mat'])
    % here we need the session file: 
    % Cell Explorer format variable containing parameters about 
    % the recording (sampling rate, shank configuration, number
    % of channels):
    % .extracellular.sr: dat sampling rate
    % .extracellular.nChannels: totel channel number
    % .extracellular.spikeGroups.channels: 1xnshank cell array
    % e.g. shk1_channel =
    % session.extracellular.spikeGroups.channels{1};
    % .general.basePath: directory to find session.dat

    Badchannels = session.channelTags.Bad.channels;% from 1 on
    Groupshanks = {[1:session.extracellular.nElectrodeGroups]};% here I use all the shanks. change according to your arrangement.
    datfile = [filename, '.dat'];
    % %% if you have sleeping periods
    % sleepstatefile = [filename, '.SleepState.states.mat'];
    % try
    %     load(sleepstatefile, 'SleepState')
    %     Sleepprd = sortrows([SleepState.ints.NREMstate;SleepState.ints.REMstate]);
    %     Sleepprd = BinPeriod(Sleepprd,5);
    %     fprintf('Skip sleep state from \n %s',sleepstatefile)
    % catch
    %     warning(sprintf('Sleep state file \n %s not exist.',sleepstatefile))
    %     Sleepprd = [];
    % end
    %%
    DenoiseDatYY(datfile,session,'SSD_path',SSD_path,'rejectChannels',Badchannels,'groupshanks',Groupshanks,'silentperiod',Sleepprd,'isoverwrite',isoverwrite_dat)

end