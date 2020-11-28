addpath('Sound_Samples/Part A/')
addpath('Sound_Samples/Part B/Violin')
addpath('Chimera code and WAV files')
addpath('BEZ2018model')
addpath('SpikeTrains_SpectroTemporal-master');
%Shouldn't have to run this more than once right?

cd BEZ2018model
mexANmodel
cd ../

clear all; clc
close all

%[test_sig test_fs]= audioread('flute_A4_normal.mp3');
[test_sig test_fs]= audioread('violin_A4_phrase_forte_arco-spiccato.mp3');

%% testing


% modFreq = 20;
% 
% 
% StimParams.fs= 100e3;
% StimParams.fm= modFreq;
% StimParams.modDepth= 1;
% StimParams.dur= 1;
% StimParams.phi_m= [];
% StimParams.dBreThresh= 40;
% StimParams.rlf_dur= .4;
% StimParams.dBSPL= 30;
% StimParams.DrivenRateTarget= 130;
% 
% ANparams.spont = 70;   % spontaneous firing rate
% ANparams.tabs= 0.6e-3; % Absolute refractory period
% ANparams.trel= 0.6e-3; % Baseline mean relative refractory period
% ANparams.cohc= 1.0;    % normal ohc function
% ANparams.cihc= 1.0;    % normal ihc function
% ANparams.species= 1;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
% ANparams.noiseType= 0;  % 1 for variable fGn; 0 for fixed (frozen) fGn
% ANparams.implnt= 0;     % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
% ANparams.dt= 1/StimParams.fs; %  time step
% 
% 
%     curCF_Hz= 1e3;
%     
%     
%     [curSAM_pos, ~]= helper.create_SAM(curCF_Hz, StimParams.fm, StimParams.fs, StimParams.modDepth, StimParams.dur, [], StimParams.phi_m);
%     
%     %     thresh_dBSPL= get_thresh_curCF(curSAM_pos(1:round(StimParams.rlf_dur*StimParams.fs)), curCF_Hz, ANparams);
%     %     curSAM_pos= gen_rescale(curSAM_pos, thresh_dBSPL + StimParams.dBreThresh);
%     
%     oa_dBSPL= helper.get_dBSPL_from_rlf(curSAM_pos(1:round(StimParams.rlf_dur*StimParams.fs)), curCF_Hz, ANparams, StimParams.DrivenRateTarget*StimParams.rlf_dur);
%     %     oa_dBSPL= 30;
% 
%     curSAM_pos= helper.gen_rescale(curSAM_pos, oa_dBSPL);
%     curSAM_neg= -curSAM_pos;
% 
% test_sig = curSAM_pos';
% test_fs = StimParams.fs;

%% Following the general structure of testANmodel_BEZ2018.m:

%TODO: Read more on these params
% model parameters
CF    = 440;   % CF in Hz; %gonna need to change this
spont = 70;   % spontaneous firing rate %SATYA CHANGED TO 70/s
tabs   = 0.6e-3; % Absolute refractory period
trel   = 0.6e-3; % Baseline mean relative refractory period
cohc  = 1.0;    % normal ohc function
cihc  = 1.0;    % normal ihc function
species = 3;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning) %read up on this tuning
noiseType = 0;  % 1 for variable fGn; 0 for fixed (frozen) fGn
implnt = 0;     % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

stimdb = 30; % stimulus intensity in dB SPL
F0 = CF; % stimulus frequency in Hz
Fmod = 40;
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 1;  % stimulus duration in seconds
rt = 2.5e-3; % rise/fall time in seconds
%ondelay = 10e-3;

% test_fs = Fs;
% t = 0:(1/test_fs):T;
% test_sig = (sin(2*pi*Fmod*t)+1).*sin(2*pi*F0*t);
test_sig = helper.gen_rescale(test_sig, stimdb);

%Modify to work with alternating polarities

nrep = 25; % number of stimulus repetitions (e.g., 50); 
psthbinwidth = 1e-4; % binwidth in seconds;
psthbins = round(psthbinwidth*Fs);  % number of psth bins per psth bin

pin = resample(test_sig, Fs, test_fs)'; 
dt=1/Fs; %  time step

% t = 0:1/Fs:T-1/Fs; % time vector
% mxpts = length(t);
% irpts = rt*Fs;
% onbin = round(ondelay*Fs);
% 
% pin = zeros(1,onbin+mxpts);
% 
% pin(onbin+1:onbin+mxpts) = sqrt(2)*20e-6*10^(stimdb/20)*sin(2*pi*F0*t); % unramped stimulus
% pin(onbin+1:onbin+irpts)= pin(onbin+1:onbin+irpts).*(0:(irpts-1))/irpts;
% pin(onbin+(mxpts-irpts):onbin+mxpts)=pin(onbin+(mxpts-irpts):onbin+mxpts).*(irpts:-1:0)/irpts;
% 

vihc_pos = model_IHC_BEZ2018(pin,CF,nrep,dt,4*T,cohc,cihc,species);
vihc_neg = model_IHC_BEZ2018(-pin,CF,nrep,dt,4*T,cohc,cihc,species);

[psth_pos, ~, ~, ~, ~,~] = model_Synapse_BEZ2018(vihc_pos,CF,nrep,dt,noiseType,implnt,spont,tabs,trel);
[psth_neg, ~, ~, ~, ~,~] = model_Synapse_BEZ2018(vihc_neg,CF,nrep,dt,noiseType,implnt,spont,tabs,trel);
% 
Psth_pos = sum(reshape(psth_pos,psthbins,length(psth_pos)/psthbins)); %
Psth_neg = sum(reshape(psth_neg,psthbins,length(psth_neg)/psthbins)); %
% binEdges= 0:1/Fs:T;
% Psth_pos = histcounts(psth_pos, binEdges);
% Psth_neg = histcounts(psth_neg, binEdges);


%TODO: 
% - Get PSD representation of both TFS and ENV 
% - Do coherence between this spectra and the spectrum of the stimulus
% (mscohere)
% 

% polarity tolerant component (ENV)
fs2 = Fs/(length(psth_pos)/length(Psth_pos));
s = (Psth_pos + Psth_neg)/2;
s = s(1:(fs2*length(test_sig)/test_fs));
%[s_psd, freqPSTH] = periodogram(s,hamming(length(s)),2048,fs2,'power'); 
[s_psd, freqPSTH]= helper.plot_dpss_psd(s,fs2); 

%FIX THIS, necessary?? Need to adjust filter params
%filtObj= helper.get_filter_fdesign('bp', [max(.1, CF-5*20) CF+5*20], fs2, 2); % second order for now


% polarity sensitve component (TFS)
d = (Psth_pos - Psth_neg)/2;
d = d(1:(fs2*length(test_sig)/test_fs)); %truncate
%d= filter(filtObj, d);

phi = sqrt(2)*rms(d)*(d./abs(hilbert(d)));
%[phi2_psd, freqs] = periodogram(phi,hamming(length(d)),2048,fs2,'centered'); 
[phi_psd, freqPSTH]= helper.plot_dpss_psd(phi,fs2); 

%% plot
simtime = length(psth_pos)/Fs;
tvect = 0:psthbinwidth:simtime-psthbinwidth;

tt= 0:1/Fs:(length(psth_pos)-1)/Fs;

px = zeros(size(psth_pos));
px(1:length(pin)) = pin;

subplot(2,1,1);
hold on
plot(tvect*1e3, Psth_pos/nrep/psthbinwidth) % Plot of estimated mean spike rate
plot(tvect*1e3, -Psth_neg/nrep/psthbinwidth) % Plot of estimated mean spike rate
hold off

ylabel('Firing Rate (/s)')
xlabel('Time (ms)')
xlim(ceil(tt([1 end])*1e3))
title('apPSTH')

subplot(2,1,2);
plot(tt,px)
ylabel('Pressure (Pa)')
xlabel('Time (ms)')
title('AN Reponse to A4 - Violin | CF = 440 Hz')


figure;
ax(1) = subplot(2,1,1);
semilogx(freqPSTH, s_psd,'b');
title('Modulation PSD')
ylabel('PSD (dB/Hz)');
xlim([0,5000])


ax(2) = subplot(2,1,2);
semilogx(freqPSTH, phi_psd,'r');
title('Carrier Freq PSD')
xlabel('Freq (Hz)');
ylabel('PSD (dB/Hz)');
xlim([0,5000])

linkaxes(ax,'x');

