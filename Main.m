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

cd Functions
mex gammatone_c.c
cd ../

%% Model Parameter Initialization:

%modelParams.CF = 0;
F0 = 440;
CF = [125, F0, 2*F0, 4e3];
dB_loss = [40, 40, 45, 50];
NFFT = 4000;
NW = 3;

%Normal Model Params
modelParams.spont = 70;
modelParams.tabs = 0.6e-3;
modelParams.trel = 0.6e-3;
modelParams.cohc = [1.0, 1.0, 1.0, 1.0]; %healthy
modelParams.cihc = [1.0, 1.0, 1.0, 1.0]; %healthy
modelParams.species = 2; % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning) %read up on this tuning
modelParams.noiseType = 0; % 1 for variable fGn; 0 for fixed (frozen) fGn
modelParams.implnt = 0; % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
modelParams.stimdb = 65;
modelParams.dur = 1; % stimulus duration in seconds, adjusted automatically to stimulus length
modelParams.reps = 25; 
modelParams.Fs = 100e3;
modelParams.psthbinwidth = 1e-4;
modelParams.buffer = 2;

%Impaired Params (just changing cohc/cihc):
[cohc_impaired,cihc_impaired,~] = fitaudiogram2(CF,dB_loss,modelParams.species);


%% Stimuli Initialization:

%Instrument
instruments = ["banjo","clarinet","flute","trombone","violin"];
%instruments = ["banjo"];
pitch = 'A4';
cond = 'normal';


%% Spectral Analysis:
l_instr = length(instruments);

bankedSig = cell(1,l_instr); 
psth_pos = cell(1,l_instr);
psth_neg = cell(1,l_instr);
env_cohere = cell(1,l_instr);
tfs_cohere = cell(1,l_instr);
s = cell(1,l_instr);
phi = cell(1,l_instr);
input_env = cell(1,l_instr);
input_tfs = cell(1,l_instr);
env_csd = cell(1,l_instr);
tfs_csd = cell(1,l_instr);
env_psd = cell(1,l_instr);
tfs_psd = cell(1,l_instr);
input_env_psd = cell(1,l_instr);
input_tfs_psd = cell(1,l_instr);

i_psth_pos = cell(1,l_instr);
i_psth_neg = cell(1,l_instr);
i_env_cohere = cell(1,l_instr);
i_tfs_cohere = cell(1,l_instr);
i_s = cell(1,l_instr);
i_phi = cell(1,l_instr);
i_env_csd = cell(1,l_instr);
i_tfs_csd = cell(1,l_instr);
i_env_psd = cell(1,l_instr);
i_tfs_psd = cell(1,l_instr);

wb = waitbar(0,'Starting Data Processing...');
for i = 1:l_instr
   
    filename = strcat(instruments(i),'_',pitch,'_',cond,'.mp3');
    
    %test
    %filename = 'SAM_test.wav';
    waitbar((i-1)/l_instr,wb,strcat('Processing- ', instruments(i),' Normal Hearing'));
    [input, input_fs] = audioread(filename);
    modelParams.dur = length(input)/input_fs;
    
    
    %Normal Hearing
    modelParams.cohc = [1.0, 1.0, 1.0, 1.0]; %healthy
    modelParams.cihc = [1.0, 1.0, 1.0, 1.0]; %healthy
    
    [bankedSig{i}] = cochlearFilterBank(input, input_fs, CF, 0); %1 for hwave rect
    [psth_pos{i}, psth_neg{i}, psth_fs] = getAP_PSTH(input, input_fs, modelParams, CF);
    [env_cohere{i}, tfs_cohere{i}, freq_cohere, s{i}, phi{i}, input_env{i}, input_tfs{i}] = getCoherence(bankedSig{i}, input_fs, psth_pos{i}, psth_neg{i}, psth_fs, NFFT, CF);
    [env_csd{i}, tfs_csd{i}, env_psd{i}, tfs_psd{i}, input_env_psd{i}, input_tfs_psd{i}, freq_SD] = getCSD(input_env{i}, input_tfs{i}, s{i}, phi{i}, NW, NFFT, psth_fs,"pmtm"); 
    
    waitbar((i-.5)/l_instr,wb,strcat('Processing- ',instruments(i),' Impaired Hearing'));

    %Impaired
    modelParams.cohc = cohc_impaired; 
    modelParams.cihc = cihc_impaired; 
    
    [i_psth_pos{i}, i_psth_neg{i}, ~] = getAP_PSTH(input, input_fs, modelParams, CF);
    [i_env_cohere{i}, i_tfs_cohere{i}, ~, i_s{i}, i_phi{i}, ~, ~] = getCoherence(bankedSig{i}, input_fs, i_psth_pos{i}, i_psth_neg{i}, psth_fs, NFFT, CF);
    [i_env_csd{i}, i_tfs_csd{i}, i_env_psd{i}, i_tfs_psd{i}, ~, ~, ~] = getCSD(input_env{i}, input_tfs{i}, i_s{i}, i_phi{i}, NW, NFFT, psth_fs,"pmtm"); 
     
end
toc

waitbar(1,wb,'Done!');
pause(0.2);
close(wb);

%% Plot Params:

instrum = 1;
CF_ind = 2;

%% apPSTH |normal/impaired| Stimulus Plot

%apPSTH
simtime = modelParams.buffer*floor(modelParams.dur);
tvect = 0:modelParams.psthbinwidth:simtime-modelParams.psthbinwidth;

tt= (0:1:(simtime*input_fs-1))/input_fs;

px = zeros(1,simtime*input_fs);
px(1:length(bankedSig{instrum}(CF_ind,:))) = bankedSig{instrum}(CF_ind,:);

figure;
subplot(3,1,1);
hold on
plot(tvect*1e3, psth_pos{instrum}(CF_ind,:)/modelParams.reps/modelParams.psthbinwidth) % Plot of estimated mean spike rate
plot(tvect*1e3, -psth_neg{instrum}(CF_ind,:)/modelParams.reps/modelParams.psthbinwidth) % Plot of estimated mean spike rate
hold off
text(75,2500,strcat('CF = ', num2str(CF(CF_ind))),'FontSize',15)

ylabel('Firing Rate (/s)')
xlim(ceil(tt([1 end])*1e3))
ylim([-2000,2000]);
title('apPSTH | Normal Hearing')
legend('(+)','( - )')
set(gca, 'FontSize',10);


subplot(3,1,2);
hold on
plot(tvect*1e3, i_psth_pos{instrum}(CF_ind,:)/modelParams.reps/modelParams.psthbinwidth) % Plot of estimated mean spike rate
plot(tvect*1e3, -i_psth_neg{instrum}(CF_ind,:)/modelParams.reps/modelParams.psthbinwidth) % Plot of estimated mean spike rate
hold off

ylabel('Firing Rate (/s)')
xlim(ceil(tt([1 end])*1e3))
ylim([-2000,2000]);
title('apPSTH | Impaired Hearing')
set(gca, 'FontSize',10);

subplot(3,1,3);
plot(tt*1e3,px,'k')
ylabel('Pressure (Pa)')
xlabel('Time (ms)')
title(strcat(instruments(instrum),' (Filtered at CF)'))
set(gca, 'FontSize',10);

set(gcf,'Position',[1200, 500, 800, 500]);
%% 