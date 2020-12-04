%Andrew Sivaprakasam

tic
%% Clearing and Adding Paths
clear all, close all

addpath('Sound_Samples/Part A/')
addpath('Sound_Samples/Part B/Tambourine')
addpath('Functions')
addpath('Sound_Samples/Part B/Violin')
addpath('Chimera code and WAV files')
addpath('BEZ2018model')
addpath('SpikeTrains_SpectroTemporal-master');

%% Compiling C Code

cd BEZ2018model
mexANmodel
cd ../ 

mex gammatone_c.c

%% Model Parameter Initialization:

%modelParams.CF = 0;
CF = [125, 440, 880, 4e3];

modelParams.spont = 70;
modelParams.tabs = 0.6e-3;
modelParams.trel = 0.6e-3;
modelParams.cohc = 1.0; %healthy
modelParams.cihc = 1.0; %healthy
modelParams.species = 3; % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning) %read up on this tuning
modelParams.noiseType = 0; % 1 for variable fGn; 0 for fixed (frozen) fGn
modelParams.implnt = 0; % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
modelParams.stimdb = 65;
modelParams.dur = 1; % stimulus duration in seconds
modelParams.reps = 25; 
modelParams.Fs = 100e3;
modelParams.psthbinwidth = 1e-4;

%% Stimuli Initialization:

%Instrument
%instruments = ["banjo","clarinet","flute","trombone","violin"];
instruments = ["violin"];
pitch = 'A4';
cond = 'normal';


%% Spectral Analysis:
l_instr = length(instruments);

bankedSig = cell(1,l_instr);
psth_pos = cell(1,l_instr);
psth_neg = cell(1,l_instr);

for i = 1:l_instr
   
    filename = strcat(instruments(i),'_',pitch,'_',cond,'.mp3');
    [input, input_fs] = audioread(filename);
    
    %Each row is a new CF
    [bankedSig{i}] = cochlearFilterBank(input, input_fs, CF, 0); %1 for hwave rect
    [psth_pos{i}, psth_neg{i}] = getAP_PSTH(input,input_fs,modelParams,CF);
    
end

%% Plotting

%apPSTH
simtime = 4*modelParams.dur;
tvect = 0:modelParams.psthbinwidth:simtime-modelParams.psthbinwidth;

hold on
plot(tvect*1e3, psth_pos{1}(3,:)/modelParams.reps/modelParams.psthbinwidth) % Plot of estimated mean spike rate
plot(tvect*1e3, -psth_neg{1}(3,:)/modelParams.reps/modelParams.psthbinwidth) % Plot of estimated mean spike rate
hold off


toc

