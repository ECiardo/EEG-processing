%% Resting State EEG Data Analysis
% This script processes resting state EEG data to analyze spectral characteristics and identify artifacts. 
% It utilizes the FieldTrip toolbox extensively.
%
% Author: Emanuele Ciardo, BCBL, March 2024
% Curator: Manuela Ruzzoli, BCBL, April 2024

%% Features:
% The script adapts the analysis described in 
% 'Alterations in rhythmic and non-rhythmic resting-state EEG activity and their link to cognition in older age' 
% (Cesnaite et al., 2023) https://pubmed.ncbi.nlm.nih.gov/36587708/
% 
% - Data Preprocessing: Filters (band-pass and notch), re-referencing, and artifact handling.
% - Independent Component Analysis (ICA): Used for artifact removal, particularly eye movements and muscle artifacts.
% - Spectral Analysis: Uses the FOOOF algorithm to model the power spectrum and separate aperiodic and oscillatory components.
% - Visualization: Generates plots of the raw and processed EEG, power spectrum, and FOOOF results.
% - Data Saving: Outputs the cleaned data and analysis results into a MATLAB file.
%
%% Usage:
% 1. Adjust the 'system' variable to match your operating system ('windows' or 'linux').
% 2. Set the appropriate paths for FieldTrip, SASICA, and your EEG data directories.
% 3. Run the script section by section to observe each phase of the processing and make necessary adjustments.
%
%% Requirements:
% - MATLAB with FieldTrip toolbox installed.
% - SASICA for MATLAB (for ICA-based artifact rejection).
% - EEG data files in FieldTrip-compatible format (e.g., .eeg, .vhdr, .vmrk).
%
%% Output:
% - A .mat file containing the cleaned EEG data, spectral analysis results, and configuration settings used.
%
%% Note:
% - Ensure that all paths are correctly set before running the script to avoid errors.
% - The script can be modified to handle different datasets or analysis parameters based on specific research needs.
%
%% Example File Names:
% - 'example_EC.eeg' for eyes-closed resting state data.
% - 'example_EO.eeg' for eyes-open resting state data.
%% -----------------------------------------------------------------------------------------------------------------------

%% General Setup
% Clear all variables, close all figures, and clear the command window for a fresh start.
clear all; close all; clc;

% Depending on the operating system, set different paths for FieldTrip and other dependencies.
system = 'windows'; % Can be 'linux' or 'windows'
switch system
    case 'windows'
        ft_path             = 'G:\MultiCRCI_DataAnalysis\fieldtrip-20230118\';  % Path to your FieldTrip folder
        SASICA_path         = 'G:\MultiCRCI_DataAnalysis\SASICA\SASICA';        % Path to your SASICA folder
        Setting.output_path = 'G:\MultiCRCI_DataAnalysis\main\output\RestingStateOutput2\'; % Output directory
        inPath              = 'H:\EEG\MULTICRCI_EEG\DATA\EEG\';                 % Input directory where EEG data is stored
        main_path           = 'G:\MultiCRCI_DataAnalysis\main\';                % Path to main project folder
    case 'linux'
        ft_path             = '/bcbl/home/public/MultiCRCI_DataAnalysis/fieldtrip-20230118';
        SASICA_path         = '/bcbl/home/public/MultiCRCI_DataAnalysis/SASICA';
        Setting.output_path = '/bcbl/home/public/MultiCRCI_DataAnalysis/main/RestingStateOutput';
        inPath              = '/bcbl/data/EEG/MULTICRCI_EEG/DATA/EEG';
        main_path           = '/bcbl/home/public/MultiCRCI_DataAnalysis/main';
end

% Add FieldTrip and SASICA directories to the MATLAB path.
addpath(ft_path);
addpath(SASICA_path);
ft_defaults;  % Initialize FieldTrip defaults
addpath(genpath(main_path));  % Add all subfolders of the main project directory to the path

% User interface to select an EEG file to process.
% [Setting.filename, Setting.input] = uigetfile({'*EC.eeg', 'EEG Files (*.EC.eeg)'}, 'Select a file', inPath); % Only EyesClosed
[Setting.filename, Setting.input] = uigetfile({'*EC.eeg;*EO.eeg', 'EEG Files (*.EC.eeg, *.EO.eeg)'}, 'Select a file', inPath); % EC and EO

% Warn the user if the file has been processed already, and ask if they
% want instead to process the first unprocessed file from the same group
if exist([Setting.output_path,Setting.filename(1:end-4),'_Processed.mat'])
    warning('The selected file has been already processed')
    go = input('Do you wish to continue with this file (0), or with the first you haven t processed (1)?');
    if go == 1
        files = dir(fullfile([Setting.input, '*',Setting.filename(end-5:end)]));
        found = false;
        for i = 1: length(files)
            focus = files(i).name;
            if ~exist([Setting.output_path,focus(1:end-4),'_Processed.mat'])
                warning(['Processing file ',focus])
                Setting.filename = focus;
                found = true;
                break
            end
        end
        if ~found
            warning('You might have processed all the files of this group. Congrats!')
        end
    end
end

% Construct full paths to the EEG data and header files.
EEG.record = append(fullfile(Setting.input), Setting.filename(1:end-4), '.eeg');
EEG.header = append(fullfile(Setting.input), Setting.filename(1:end-4), '.vhdr');
EEG.event  = append(fullfile(Setting.input), Setting.filename(1:end-4), '.vmrk');

% Read the EEG header to get information about the data.
cfg = [];
cfg.dataset    = [EEG.record];
cfg.headerfile = [EEG.header];
cfg.eventfile  = [EEG.event];
    hdr = ft_read_header(cfg.headerfile); % variable not used

% Read sensor locations necessary for later analysis steps.
elec = ft_read_sens('G:\MultiCRCI_DataAnalysis\fieldtrip-20230118\template\electrode\easycap-M1.txt');

%% Preprocessing
% Set preprocessing options, including reference electrodes and filtering.
cfg.reref      = 'yes';
cfg.refchannel = {'M1', 'M2'};      % Reference to mastoid channels
cfg.channel    = {'all','-Audio'};  % Select all but audio channel

% Configure band-pass filtering from 1-45 Hz using a 4th order Butterworth filter.
cfg.padding   = 10; %add a padding before filtering
cfg.bpfilter  = 'yes';
cfg.bpfreq    = [1 45]; % same as Cesnaite et at. 2023
cfg.bpfiltord = 4;

% Configure notch filtering to remove line noise at 50 Hz.
cfg.dftfilter = 'yes';
cfg.dftfreq   = 50;

% Apply preprocessing to the data.
    data1 = ft_preprocessing(cfg);

% Assign electrode information back to the data structure.
data1.elec = elec; % this can be useful for later analysis. I added it as I was trying to make work a function I eventually discarded
Setting.elect = elec; %  % save info

% Discard from the data the first few seconds to account for initial movement or setup.
discardStart = 10; % Time to exclude at the beginning, in seconds
cfg = [];
cfg.latency = [discardStart floor((max(data1.time{1})))]; % up to the lower integer in the data
    data2 = ft_selectdata(cfg, data1);
Setting.discardStart = discardStart; % save info

% Epoch the data into 1-second segments without overlap. This helps in handling artefacts effectively.
cfg = [];
cfg.length  = 1;
cfg.overlap = 0;
    data3 = ft_redefinetrial(cfg, data2);

%% Manual Artifact Rejection
% Browse the data manually to identify and reject artifacts.
cfg = [];
cfg.blocksize  = 8;
cfg.continuous = 'yes';
cfg.channel    = {'all'}; % Leave all channels to be visually aided in detecting non EEG activity
    cfg = ft_databrowser(cfg, data3);
ManualArt = cfg.artfctdef.visual.artifact;
Setting.ManualArt = ManualArt;  % save info

% Reject the epochs containing artifacts
cfg.artfctdef.reject = 'complete';
    data4 = ft_rejectartifact(cfg, data3);

%% Artifacts by rejectvisual: select bad trials or channels from the summary plot, BE CONSERVATIVE!
%Just in case there is a bad channel you need to remove
% Write down the excluded trials and channels
cfg = [];
cfg.layout = 'acticap-64ch-standard2.mat';
cfg.method = 'summary'; %'channel'; %'summary'; 'trial';
% cfg.channel = 'EEG';
cfg.keepchannel = 'no';
    data4 = ft_rejectvisual(cfg, data4);

%% (Manual) ICA artifact correction
% Perform Independent Component Analysis (ICA) to iden1tify components related to eye movements or other artifacts.
% in order to get rid of blink artefacts and external noise, the
% independent component analysis can be useful. But Beware, don't exclude
% to much! Whereas the ICA can split the signal into independent
% components, this does not mean that one component does not contain
% information from more than one source! 
% Have a look at ICA Lecture in http://mikexcohen.com/lectures.html
clearvars data1 data2 data3
tic; % Start timer for performance measurement
cfg = [];
%EOG included in ICA, see https://www.researchgate.net/post/Include_EOG_in_EEG_ICA_or_not
cfg.channel = {'all', '-M1', '-M2', '-Audio'};  % Include all channels
% cfg.channel = {'all', '-M1', '-M2', '-RVEOG','-LHEOG','-RHEOG','-Audio'};

cfg.method = 'runica';  % Use the FastICA algorithm
cfg.randomseed = [1516]; % to make sure that component N will be the same the next time you run it
cfg.numcomponent = 30; % Limit the number of components. 
% the results with all ICs
    comp = ft_componentanalysis(cfg, data4);
toc % Stop timer and print elapsed time

%% Take a look at the components, play around with the way you want to see the components  
% % % FFT of the components
% for cc = 1:20
% CompNumber = cc; %write the number of the component you want the FFT 
% 
% cfg = [];
% cfg.channel = comp.label{CompNumber};
% cfg.method = 'mtmfft';
% cfg.pad = 'nextpow2';
% cfg.output      = 'pow';
% cfg.taper      = 'hanning';
% cfg.foilim     = [1 50];
%     freq = ft_freqanalysis(cfg, comp);
% subplot (5,4,cc); plot(freq.freq, freq.powspctrm); title(['FFT Comp ', num2str(cc)]); hold on 
% end
% 
% % %% Average components
% % Topoplot
% cfg           = [];
% cfg.colormap  = 'jet';
% cfg.component = 1:20;
% cfg.comment   = 'no';
% cfg.layout    = 'acticap-64ch-standard2.mat'; %'biosemi64.lay'; % If you use a function that requires plotting of topographical information you need to supply the function with the location of your channels
% figure;     
%     ft_topoplotIC(cfg, comp);
%     
% % Single component, time-course
% % cfg = [];
% % cfg.layout    = 'biosemi64.lay';
% % cfg.viewmode = 'component';
% %     ft_databrowser(cfg, comp);    
% 
% % Calculate the average comp
% % Doesn't make a lot of sense because it is not an event-related design
% % cfg = [];
% % comp_avg = ft_timelockanalysis(cfg, comp);  
% % % Single component, grandaverage time-course
% % cfg = [];
% % cfg.layout    = 'acticap-64ch-standard2.mat'; %'biosemi64.lay';
% % cfg.viewmode = 'component';
% %     ft_databrowser(cfg, comp_avg);
% 
% %% ----------- ATTENTION HERE, MUST WRITE THE COMPONENTS to EXCLUDE -------
% %% And take out the ones that are clearly artefacts
% close all;
% cfg = [];
% cfg.component = [5,8]; % I focus on blinks, muscle artefacts, envirormental noise, heart
%     data5 = ft_rejectcomponent(cfg,comp, data4);
% Setting.RemovedICA = cfg.component; % save info
% %% ------------------------------------------------------------------------

%% SASICA
% Using SASICA for a computer-assisted selection of bad ICs
% Please read this paper (https://pubmed.ncbi.nlm.nih.gov/25791012/) from
% the developers of SASICA which explain how to identify bad components.
% Note that the user will receive suggestions from both SASICA and Mohamad
% script about potential bad ICs, but still will have to look at them and
% make their own decision

% create a data file for SASICA (it needs to have the same channels as the
% ones specified while performing ICA
cfg = [];
cfg.channel = {'all', '-M1', '-M2', '-Audio'};
dataSASICA  = ft_selectdata(cfg,data4);

% set up SASICA params
cfg = []; 
    cfg = ft_SASICA('getdefs'); % Initialize cfg with defaults values

% Specify which analysis to run in SASICA; I allowed all them all
cfg.layout                = 'easycapM11';
cfg.autocorr.enable       = 1;
cfg.focalcomp.enable      = 1;
cfg.SNR.enable            = 1;
cfg.EOGcorr.enable        = 1;
cfg.EOGcorr.Veogchannames = {'RVEOG'};
cfg.EOGcorr.Heogchannames = {'RHEOG', 'LHEOG'};
cfg.chancorr.enable       = 1;
cfg.FASTER.enable         = 1;
cfg.FASTER.blinkchanname  = {'RVEOG'};
% cfg.ADJUST.enable         = 1;


close all
    cfg_ICA = ft_SASICA(cfg, comp, dataSASICA);

% Apply ICA component rejection to clean the data.
cfg = [];

%%%% problem here
cfg.component = double(cfg_ICA.reject.gcompreject);
    data5 = ft_rejectcomponent(cfg, comp, data4);

% %% Visualization of the results pre- and post-ICA cleaning
% 
% 
% 
% close all;
% figure('units', 'normalized', 'OuterPosition', [0, 0.2, 1, 0.6]);
% hold on;
% 
% % Concatenate data for visual comparison
% 
% s1 = subplot(211);
% hold on;
% chan = 16;  % Example frontal channel number
% showData = [];
% for i = 1: size(data4.trial,2)
%     showData = [showData, data4.trial{i}(chan,:)];
% end
% plot([1:length(showData)]/data4.fsample, showData, 'DisplayName', 'Original');
% 
% showData = [];
% for i = 1: size(data5.trial,2)
%     showData = [showData, data5.trial{i}(chan,:)];
% end
% plot([1:length(showData)]/data5.fsample, showData, 'DisplayName', 'ICA-cleaned');
% xlabel('Time (sec)');
% ylabel('Power');
% legend show;
% title('Frontal Channel')
% xlim([100,110]);  % Limit x-axis for better visibility
% 
% s2 = subplot(212);
% hold on;
% chan = 59;  % Example occipital channel number
% showData = [];
% for i = 1: size(data4.trial,2)
%     showData = [showData, data4.trial{i}(chan,:)];
% end
% plot([1:length(showData)]/data4.fsample, showData, 'DisplayName', 'Original');
% 
% showData = [];
% for i = 1: size(data5.trial,2)
%     showData = [showData, data5.trial{i}(chan,:)];
% end
% plot([1:length(showData)]/data5.fsample, showData, 'DisplayName', 'ICA-cleaned');
% xlabel('Time (sec)');
% ylabel('Power');
% legend show;
% title('Occipital Channel')
% xlim([100,110]);  % Limit x-axis for better visibility
% linkaxes([s1,s2])

%% Artifacts by rejectvisual: select bad trials or channels from the summary plot, BE CONSERVATIVE!
% Write down the excluded trials and channels
cfg = [];
cfg.layout = 'acticap-64ch-standard2.mat';
cfg.method = 'summary'; %'channel'; %'summary'; 'trial';
cfg.channel = 'EEG';
cfg.keepchannel = 'no';
    CleanData = ft_rejectvisual(cfg, data5);

%% Clear same variables
clearvars -except data4 data5 CleanData Setting ft_path inPath main_path 

%% Spectral analysis using FOOOF (fitting oscillations and one over f)
% Configure spectral analysis parameters
cfg = [];
cfg.foi       = [1:0.1:40];  % Range of frequencies of interest
cfg.channel   = {'eeg', '-M1', '-M2','-RVEOG', '-LHEOG', '-RHEOG'}; 
cfg.taper     = 'hanning';
cfg.pad       = 4;  % Padding for spectral analysis
cfg.tapsmofrq = 2;
cfg.method    = 'mtmfft';
cfg.output    = 'fooof_aperiodic';  % Output parameter set to 'fooof_aperiodic'
cfg.fooof.aperiodic_mode = 'fixed'; % Options: 'fixed' or 'knee'

% Execute frequency analysis
%     fractal = ft_freqanalysis(cfg, data5);
fractal = ft_freqanalysis(cfg, CleanData);
cfg.output = 'pow';
%     original = ft_freqanalysis(cfg, data5);
original = ft_freqanalysis(cfg, CleanData);

% Subtract the fractal component from the power spectrum to isolate oscillatory components
cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'x2-x1';  % Subtract the aperiodic component
    oscillatory = ft_math(cfg, fractal, original);

% Alternatively, compute the oscillatory component as the quotient of the power spectrum and the fractal component
% The original paper calls for a subtraction method, but we could consider
% adding a remark on the division method.
cfg.operation = 'x2./x1';  % Quotient operation
    oscillatory_ratio = ft_math(cfg, fractal, original);

%% Identify and flag channels with poor FOOOF model fits
cutoff = 0.80; % cutoff parameter. 0.8 is what used in Cesnaite et al., 2023
flag = []; % vector of channels with poor fit
for i = 1:size(fractal.fooofparams, 1)
    if fractal.fooofparams(i).r_squared < cutoff
        flag = [flag, i];
    end
end
% Write on screen the bad electrodes
fractal.label{flag(1:end), 1}  

%% Visualization of FOOOF results
% Plot the original, aperiodic, and oscillatory components in a multi-panel figure
close all;
figure('units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
trsp = 0.2;  % Transparency for overlay plots

subplot(221) % Original power spectrum
hold on;
f = fractal.freq;
focus = original.powspctrm;

% select frontal channels
temp.frontal = {'Fp1','Fp2','AF3','AF4','AF7','AF8','AFz','F1','F2','F3','F4','F5','F6','F7','F8','Fz'};
frontChan = [];
for chan = temp.frontal
   frontChan = [frontChan, find(strcmp (CleanData.label,chan))];
end

temp.occipital = {'P1','P2','P3','P4','P5','P6','P7','P8','Pz','PO3','PO4','PO7','PO8','POz','O1','O2','Oz'};
occiChan = [];
for chan = temp.occipital
   occiChan = [occiChan, find(strcmp (CleanData.label,chan))];
end

plot(f, mean(focus(frontChan,:)), 'Color', [0, 0, 1], 'LineWidth',2); %% 'Fp1','Fp2','AF3','AF4','AF7','AF8','AFz','F1','F2','F3','F4','F5','F6','F7','F8','Fz'
plot(f, mean(focus(occiChan,:)), 'Color', [1, 0, 0], 'LineWidth',2); %% 'P1','P2','P3','P4','P5','P6','P7','P8','Pz','PO3','PO4','PO7','PO8','POz','O1','O2','Oz'

plot(f, focus(frontChan,:)', 'Color', [0, 0, 1, trsp]);
plot(f, focus(occiChan,:)', 'Color', [1, 0, 0, trsp]);
xlim([min(f), max(f)]);
xlabel('Frequencies (Hz)');
ylabel('Amplitude');
title('Original Spectrogram');
legend('Frontal', 'Occipital');

subplot(222) % Aperiodic component
hold on;
focus = fractal.powspctrm;
plot(f, mean(focus(frontChan,:)), 'Color', [0, 0, 1], 'LineWidth',2);
plot(f, mean(focus(occiChan,:)), 'Color', [1, 0, 0], 'LineWidth',2);
plot(f, focus(frontChan,:)', 'Color', [0, 0, 1, trsp]);
plot(f, focus(occiChan,:)', 'Color', [1, 0, 0, trsp]);
xlim([min(f), max(f)]);
xlabel('Frequencies (Hz)');
ylabel('Amplitude');
title('Aperiodic Component');

subplot(223) % Oscillatory component
hold on;
focus = oscillatory.powspctrm;
plot(f, mean(focus(frontChan,:)), 'Color', [0, 0, 1], 'LineWidth',2);
plot(f, mean(focus(occiChan,:)), 'Color', [1, 0, 0], 'LineWidth',2);
plot(f, focus(frontChan,:)', 'Color', [0, 0, 1, trsp]);
plot(f, focus(occiChan,:)', 'Color', [1, 0, 0, trsp]);
xlim([min(f), max(f)]);
xlabel('Frequencies (Hz)');
ylabel('Amplitude');
title('Oscillatory Component');

subplot(224) % Plot R-squared values for FOOOF fit across channels
hold on;
for i = 1:length(fractal.label)
    plot(i, fractal.fooofparams(i).r_squared, 'o', 'Color', 'blue');
end
plot([1, length(fractal.label)], [0.8, 0.8], '--', 'Color', 'black');
title('R squared');
xlabel('Channels');

%% Find peaks in the region of interest for the oscillatory and original data
EOI_IAF = {'P1','P2','P3','P4','P5','P6','P7','P8','Pz','PO3','PO4','PO7','PO8','POz','O1','O2','Oz'};

cfg = [];
cfg.channel = EOI_IAF;
cfg.avgoverchan = 'yes';
plot_data_oscillatory = ft_selectdata(cfg, oscillatory);
plot_data_original = ft_selectdata(cfg, original );

%% Find peaks - oscillatory
[pks_oscillatory, locs_oscillatory] = findpeaks(plot_data_oscillatory.powspctrm,plot_data_oscillatory.freq,'SortStr','descend');

% Define the alpha band (8-14)
alpha_pks_oscillatory = pks_oscillatory((locs_oscillatory<14)&(locs_oscillatory>8));
alpha_locs_oscillatory = locs_oscillatory((locs_oscillatory<14)&(locs_oscillatory>8));

if isempty(alpha_locs_oscillatory)
    IAF_oscillatory = NaN;
    IAF_ampli_oscillatory = NaN;
    fprintf('No peaks in the 8-14 Hz range detected. \n');
else
    IAF_oscillatory = alpha_locs_oscillatory(1);
    IAF_ampli_oscillatory = alpha_pks_oscillatory(1);
    fprintf('Oscillatory: Individual Alpha Frequency for this participant is: %4.2f Hz \n',IAF_oscillatory);
    fprintf('Oscillatory: The amplitude of the Individual Alpha Peak for this participant is: %4.2f \n',IAF_ampli_oscillatory);
end

    %% Find peaks - original
[pks_original, locs_original] = findpeaks(plot_data_original.powspctrm,plot_data_original.freq,'SortStr','descend');

% Define the alpha band (8-14)
alpha_pks_original = pks_original((locs_original<14)&(locs_original>8));
alpha_locs_original = locs_original((locs_original<14)&(locs_original>8));

if isempty(alpha_locs_original)
    IAF_original = NaN;
    IAF_ampli_original = NaN;
    fprintf('No peaks in the 8-14 Hz range detected. \n');
else
    IAF_original = alpha_locs_original(1);
    IAF_ampli_original = alpha_pks_original(1);
    fprintf('Original: Individual Alpha Frequency for this participant is: %4.2f Hz \n',IAF_original);
    fprintf('Original: The amplitude of the Individual Alpha Peak for this participant is: %4.2f \n',IAF_ampli_original);
end

%% Store Info IAF
clear IAF
IAF.original.pks = IAF_original;
IAF.original.locs = IAF_ampli_original;
IAF.oscillatory.pks = IAF_oscillatory;
IAF.oscillatory.locs = IAF_ampli_oscillatory;

%% Save the Results
% Define the filename for saving the results based on the original EEG file
outputFilename = fullfile(Setting.output_path, [Setting.filename(1:end-4),'_Processed']);


% Save the cleaned data
save(outputFilename, 'data5', 'fractal', 'oscillatory', 'oscillatory_ratio', 'original', 'IAF', 'Setting', 'CleanData');

%% Generate plots for topography and Individual alpha peak
figure;
cfg             = [];
cfg.layout = 'acticap-64ch-standard2.mat';
% cfg.xlim        = [8 14];
cfg.xlim        = IAF.original.pks + [-1 1];
cfg.style       = 'straight';
cfg.colorbar    = 'SouthOutside';
cfg.parameter   = 'powspctrm';
cfg.comment     = 'auto';
    ft_topoplotER(cfg,original);
title('PSD (8-14 Hz), original');
hold on
cfg.xlim        = IAF.oscillatory.pks + [-1 1];
    ft_topoplotER(cfg,oscillatory);
title('PSD (8-14 Hz), oscillatory');
%     colormap('jet');


