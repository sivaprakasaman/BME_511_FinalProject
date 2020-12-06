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
dB_stim = 65;

%Normal Model Params
modelParams.spont = 70;
modelParams.tabs = 0.6e-3;
modelParams.trel = 0.6e-3;
modelParams.cohc = [1.0, 1.0, 1.0, 1.0]; %healthy
modelParams.cihc = [1.0, 1.0, 1.0, 1.0]; %healthy
modelParams.species = 2; % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning) %read up on this tuning
modelParams.noiseType = 0; % 1 for variable fGn; 0 for fixed (frozen) fGn
modelParams.implnt = 0; % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
modelParams.stimdb = dB_stim;
modelParams.dur = 1; % stimulus duration in seconds, adjusted automatically to stimulus length
modelParams.reps = 75; 
modelParams.Fs = 100e3;
modelParams.psthbinwidth = 1e-4;
modelParams.buffer = 2;

%Impaired Params (just changing cohc/cihc):
[cohc_impaired,cihc_impaired,~] = fitaudiogram2(CF,dB_loss,modelParams.species);


%% Stimuli Initialization:

%Instrument
%instruments = ["banjo","clarinet","flute","trombone","violin"];
%instruments = ["banjo"];
instruments = ["SAM Tone"];

%ARTICULATION
% instruments = ["Spiccato","Martele"];
% articulations = ["violin_A4_phrase_forte_arco-spiccato.mp3","violin_A4_phrase_forte_arco-martele.mp3"];

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

dur = cell(1,l_instr);

wb = waitbar(0,'Starting Data Processing...');
for i = 1:l_instr
   
    %instrument comparison
    %filename = strcat(instruments(i),'_',pitch,'_',cond,'.mp3');
    
    %articulation comparison
    %filename = articulations(i);
    
    %test
    filename = 'SAM_test.wav';
    waitbar((i-1)/l_instr,wb,strcat('Processing- ', instruments(i),' Normal Hearing'));
    [input, input_fs] = audioread(filename);
    
    dur{i} = length(input)/input_fs;
    modelParams.dur = dur{i};
    
    %Normal Hearing
    modelParams.cohc = [1.0, 1.0, 1.0, 1.0]; %healthy
    modelParams.cihc = [1.0, 1.0, 1.0, 1.0]; %healthy
    modelParams.stimdb = dB_stim;

    [bankedSig{i}] = cochlearFilterBank(input, input_fs, CF, 0); %1 for hwave rect
    [psth_pos{i}, psth_neg{i}, psth_fs] = getAP_PSTH(input, input_fs, modelParams, CF);
    [env_cohere{i}, tfs_cohere{i}, freq_cohere, s{i}, phi{i}, input_env{i}, input_tfs{i}] = getCoherence(bankedSig{i}, input_fs, psth_pos{i}, psth_neg{i}, psth_fs, NFFT, CF);
    [env_csd{i}, tfs_csd{i}, env_psd{i}, tfs_psd{i}, input_env_psd{i}, input_tfs_psd{i}, freq_SD] = getCSD(input_env{i}, input_tfs{i}, s{i}, phi{i}, NW, NFFT, psth_fs,"pmtm"); 
    
    waitbar((i-.5)/l_instr,wb,strcat('Processing- ',instruments(i),' Impaired Hearing'));

    %Impaired
    modelParams.cohc = cohc_impaired; 
    modelParams.cihc = cihc_impaired; 
    modelParams.stimdb = dB_stim + mean(dB_loss); %amplified

    [i_psth_pos{i}, i_psth_neg{i}, ~] = getAP_PSTH(input, input_fs, modelParams, CF);
    [i_env_cohere{i}, i_tfs_cohere{i}, ~, i_s{i}, i_phi{i}, ~, ~] = getCoherence(bankedSig{i}, input_fs, i_psth_pos{i}, i_psth_neg{i}, psth_fs, NFFT, CF);
    [i_env_csd{i}, i_tfs_csd{i}, i_env_psd{i}, i_tfs_psd{i}, ~, ~, ~] = getCSD(input_env{i}, input_tfs{i}, i_s{i}, i_phi{i}, NW, NFFT, psth_fs,"pmtm"); 
     
end
toc

waitbar(1,wb,'Done!');
pause(0.2);
close(wb);

%% Audiogram:

audiogram = figure;

hold on
plot(CF, zeros(1,length(dB_loss)),'o-k','LineWidth',4,'MarkerSize',15);
plot(CF, dB_loss,'o-r','LineWidth',4,'MarkerSize',15);
hold off

xlim([min(CF)-20,max(CF)+100]);
ylim([-5, max(dB_loss)+5]);
set(gca, 'Ydir', 'reverse')
set(gca, 'XScale', 'log')
set(gca, 'FontSize', 12);
xticks(CF);
xticklabels(split(num2str(CF)));
xtickangle(45)
title('Audiogram');
xlabel('Frequency (Hz)');
ylabel('dB Hearing Loss');
legend('Normal','Impaired','Location','southwest');
grid on;
box on;

%% Plot Params:

instrum = 1;
CF_ind = 2;

%% apPSTH |normal/impaired| Stimulus Plot

%apPSTH
simtime = modelParams.buffer*ceil(dur{instrum});
tvect = 0:modelParams.psthbinwidth:simtime-modelParams.psthbinwidth;

tt= (0:1:(simtime*input_fs-1))/input_fs;

px = zeros(1,simtime*input_fs);
px(1:length(bankedSig{instrum}(CF_ind,:))) = bankedSig{instrum}(CF_ind,:);

tenv =  (0:1:(simtime*psth_fs-1))/psth_fs;

px_env = zeros(1,simtime*psth_fs);
px_env(1:length(input_env{instrum}(:,CF_ind))) = input_env{instrum}(:,CF_ind);

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
hold on
plot(tt*1e3, px,'k');
plot(tenv*1e3, px_env,'r','LineWidth',1.5);
hold off
ylabel('Pressure (Pa)')
xlabel('Time (ms)')
title(strcat(instruments(instrum),' (Filtered at CF)'))
legend('Stimulus','Envelope')
set(gca, 'FontSize',10);

set(gcf,'Position',[1200, 500, 800, 500]);

%% Coherence:

instrum = instrum;
CF_ind = CF_ind;

figure;
subplot(2,1,1);
hold on
plot(freq_cohere, tfs_cohere{instrum}(:,CF_ind),'LineWidth',2.5);
plot(freq_cohere, i_tfs_cohere{instrum}(:,CF_ind),'LineWidth',2.5);
hold off
%set(gca, 'XScale', 'log')
set(gca, 'FontSize',10);
xlim([CF(CF_ind)-100, CF(CF_ind)+100]);
ylim([0,1]);
title('TFS Mag. Squared Coherence')
legend('Normal','Impaired');
ylabel('Coherence');

text(CF(CF_ind)-65,.85,strcat('CF = ', num2str(CF(CF_ind))),'FontSize',15,'HorizontalAlignment','center')
text(CF(CF_ind)-65,.75,instruments(instrum),'FontSize',15,'HorizontalAlignment','center')

box on;
grid on;


subplot(2,1,2);
plot(freq_cohere,tfs_cohere{instrum}(:,CF_ind)-i_tfs_cohere{instrum}(:,CF_ind),'k','LineWidth',2.5);
title('Normal Hearing - Impaired Hearing TFS Coherence')
legend('NH minus HL Coherence');
xlim([CF(CF_ind)-100, CF(CF_ind)+100]);
ylim([-1,1])
xlabel('Frequency (Hz)');
ylabel('Difference');
set(gca, 'FontSize',10);
box on;
grid on;
set(gcf,'Position',[1200, 300, 800, 600]);

figure;
subplot(2,1,1);
hold on
plot(freq_cohere,env_cohere{instrum}(:,CF_ind),'LineWidth',2.5);
plot(freq_cohere, i_env_cohere{instrum}(:,CF_ind),'LineWidth',2.5);
hold off
set(gca, 'XScale', 'log')
set(gca, 'FontSize',10);
%xlim([CF(CF_ind)-100, CF(CF_ind)+100]);
ylim([0,1]);
title('ENV Mag. Squared Coherence')
legend('Normal','Impaired');
xlabel('Frequency (Hz)');
ylabel('Coherence');

text(5,.85,strcat('CF = ', num2str(CF(CF_ind))),'FontSize',15,'HorizontalAlignment','center')
text(5,.75,instruments(instrum),'FontSize',15,'HorizontalAlignment','center')

box on;
grid on;


subplot(2,1,2);
plot(freq_cohere,env_cohere{instrum}(:,CF_ind)-i_env_cohere{instrum}(:,CF_ind),'k','LineWidth',2.5);
legend('NH minus HL Coherence');
%xlim([CF(CF_ind)-100, CF(CF_ind)+100]);
title('Normal Hearing - Impaired Hearing ENV Coherence')
xlabel('Frequency (Hz)');
ylabel('Difference');
ylim([-1,1])
set(gca, 'XScale', 'log')
box on;
grid on;
set(gcf,'Position',[2000, 300, 800, 600]);




figure;

hold on
for i = 1:length(instruments)
    plot(freq_cohere,tfs_cohere{i}(:,CF_ind)-i_tfs_cohere{i}(:,CF_ind),'LineWidth',2.5);
end
hold off
%set(gca, 'XScale', 'log')
legend(instruments);
xlim([CF(CF_ind)-100, CF(CF_ind)+100]);
title('Normal Hearing - Impaired Hearing TFS Coherence')
xlabel('Frequency (Hz)');
ylabel('Difference');
ylim([-1,1])
box on;
grid on;
set(gcf,'Position',[1200, 2000, 800, 600]);
text(CF(CF_ind)-65,.85,strcat('CF = ', num2str(CF(CF_ind))),'FontSize',15,'HorizontalAlignment','center')


figure;
hold on
for i = 1:length(instruments)
   
    plot(freq_cohere,env_cohere{i}(:,CF_ind)-i_env_cohere{i}(:,CF_ind),'LineWidth',2.5);
    
end
hold off
set(gca, 'XScale', 'log')
title('Normal Hearing - Impaired Hearing ENV Coherence')
xlabel('Frequency (Hz)');
ylabel('Difference');
ylim([-1,1])
legend(instruments);

box on;
grid on;
text(5,.85,strcat('CF = ', num2str(CF(CF_ind))),'FontSize',15,'HorizontalAlignment','center')
set(gcf,'Position',[2000, 2000, 800, 600]);

%%
