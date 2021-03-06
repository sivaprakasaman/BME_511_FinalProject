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

[test_sig test_fs]= audioread('violin_A4_normal.mp3');
%[test_sig test_fs]= audioread('violin_A4_phrase_forte_arco-spiccato.mp3');

%% testing
% 

modFreq = 20;


StimParams.fs= 100e3;
StimParams.fm= modFreq;
StimParams.modDepth= 1;
StimParams.dur= 1;
StimParams.phi_m= [];
StimParams.dBreThresh= 40;
StimParams.rlf_dur= .4;
StimParams.dBSPL= 30;
StimParams.DrivenRateTarget= 130;

ANparams.spont = 70;   % spontaneous firing rate
ANparams.tabs= 0.6e-3; % Absolute refractory period
ANparams.trel= 0.6e-3; % Baseline mean relative refractory period
ANparams.cohc= 1.0;    % normal ohc function
ANparams.cihc= 1.0;    % normal ihc function
ANparams.species= 2;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
ANparams.noiseType= 0;  % 1 for variable fGn; 0 for fixed (frozen) fGn
ANparams.implnt= 0;     % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
ANparams.dt= 1/StimParams.fs; %  time step


    curCF_Hz= 1e3;
    
    
    [curSAM_pos, ~]= helper.create_SAM(curCF_Hz, StimParams.fm, StimParams.fs, StimParams.modDepth, StimParams.dur, [], StimParams.phi_m);
    
    %     thresh_dBSPL= get_thresh_curCF(curSAM_pos(1:round(StimParams.rlf_dur*StimParams.fs)), curCF_Hz, ANparams);
    %     curSAM_pos= gen_rescale(curSAM_pos, thresh_dBSPL + StimParams.dBreThresh);
    
    oa_dBSPL= helper.get_dBSPL_from_rlf(curSAM_pos(1:round(StimParams.rlf_dur*StimParams.fs)), curCF_Hz, ANparams, StimParams.DrivenRateTarget*StimParams.rlf_dur);
    %     oa_dBSPL= 30;

    curSAM_pos= helper.gen_rescale(curSAM_pos, oa_dBSPL);
    curSAM_neg= -curSAM_pos;

test_sig = curSAM_pos';
test_fs = StimParams.fs;

%% Generate Spikes:

%TODO: Read more on these params
% model parameters
CF    = 440;   % CF in Hz; %gonn5a need to change this
spont = 70;   % spontaneous firing rate %SATYA CHANGED TO 70/s
tabs   = 0.6e-3; % Absolute refractory period
trel   = 0.6e-3; % Baseline mean relative refractory period
cohc  = .1;    % normal ohc function
cihc  = .1;    % normal ihc function
species = 3;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning) %read up on this tuning
noiseType = 0;  % 1 for variable fGn; 0 for fixed (frozen) fGn
implnt = 0;     % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

stimdb = 65; % stimulus intensity in dB SPL
F0 = CF; % stimulus frequency in Hz
Fmod = 40;
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 1;  % stimulus duration in seconds
rt = 2.5e-3; % rise/fall time in seconds
%ondelay = 10e-3;

% test_fs = Fs;
% t = 0:(1/test_fs):T;
% test_sig = (sin(2*pi*Fmod*t)+1).*sin(2*pi*F0*t);
% test_sig = test_sig';
test_sig = helper.gen_rescale(test_sig, stimdb);

%Modify to work with alternating polarities

nrep = 25; % number of stimulus repetitions (e.g., 50); 
psthbinwidth = 1e-4; % binwidth in seconds;
psthbins = round(psthbinwidth*Fs);  % number of psth bins per psth bin

pin = resample(test_sig, Fs, test_fs)';
pin = pin(1:T*Fs);
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
%% Spectral Analyses

%TODO: 
% - Get PSD representation of both TFS and ENV 
% - Do coherence between this spectra and the spectrum of the stimulus
% (mscohere)
% 

NW = 2;
NFFT = 4e3;
lp_co = 500;
lp_ord = 6;

fs2 = Fs/(length(psth_pos)/length(Psth_pos));
input = resample(pin,fs2,Fs);

%Hilbert Env/filtering:
%hilb_env = abs(hilbert(input));

% [B,A] = butter(lp_ord,lp_co/(Fs/2));
% hilb_env = filter(B,A,hilb_env);
hilb_env = envelope(input,75,'peak');
[sig_psd_pmtm, freq_pmtm_sig] = pmtm(input,NW,NFFT,fs2);
[sig_psd_env, freq_pmtm_env] = pmtm(hilb_env,NW,NFFT,fs2);
[sig_psd_tfs, freq_pmtm_tfs] = pmtm(cos(angle(hilbert(input))),NW,NFFT,fs2);

% polarity tolerant component (ENV)
s = (Psth_pos + Psth_neg)/2;
s = s(1:(fs2*T));
%[s_psd, freqPSTH] = periodogram(s,hamming(length(s)),2048,fs2,'power'); 
[s_psd, freqPSTH]= helper.plot_dpss_psd(s,fs2); 
[s_psd_w, freq_w] = pwelch(s, [], [] ,NFFT, fs2);
%try PMTM
[s_psd_pmtm, freq_pmtm] = pmtm(s,NW, NFFT,fs2);

%FIX THIS, necessary?? Need to adjust filter params
% filtObj= helper.get_filter_fdesign('bp', [max(.1, CF-5*20) CF+5*20], fs2, 2); % second order for now


% polarity sensitve component (TFS)
d = (Psth_pos - Psth_neg)/2;
d = d(1:(fs2*T)); %truncate
%d= filter(filtObj, d);

phi = sqrt(2)*rms(d)*(d./abs(hilbert(d)));
%[phi2_psd, freqs] = periodogram(phi,hamming(length(d)),2048,fs2,'centered'); 
[phi_psd, freqPSTH]= helper.plot_dpss_psd(phi,fs2); 
[phi_psd_w, freq_w] = pwelch(phi, [], [] ,4e3, fs2);
 
[phi_psd_pmtm, freq_pmtm] = pmtm(phi,NW, NFFT,fs2);




sig_psd_pmtm = 10*log10(sig_psd_pmtm);
sig_psd_env = 10*log10(sig_psd_env);
sig_psd_tfs = 10*log10(sig_psd_tfs);

s_psd_w = 10*log10(s_psd_w);
s_psd_pmtm = 10*log10(s_psd_pmtm);
phi_psd_w = 10*log10(phi_psd_w);
phi_psd_pmtm = 10*log10(phi_psd_pmtm);

%% Coherence & Cross-Spectral Density (Depends on phase, so cannot use PSD outputs)

[TFS_coherence f_coherence] = mscohere(input,phi, [],[],NFFT,fs2);
[ENV_coherence f_coherence] = mscohere(hilb_env,s,[],[],NFFT,fs2);

% [TFS_coherence f_coherence] = cpsd(input,phi, [],[],NFFT,fs2);
% [ENV_coherence f_coherence] = cpsd(s,hilb_env, [],[],NFFT,fs2);

%this might be "weird". Figure out a workaround

WINDOW = hamming(length(input)/3);
OVERLAP = .75;

tfs_corr_hilb_tfs = xcorr(input,phi);
%P_tfs_hilb_tfs = pmtm(tfs_corr_hilb_tfs,NW, NFFT,fs2);
P_tfs_hilb_tfs = pwelch(tfs_corr_hilb_tfs,WINDOW,OVERLAP,NFFT, fs2);
P_tfs_hilb_tfs = 10*log10(P_tfs_hilb_tfs);

env_corr_hilb_env = xcorr(hilb_env,s);
%P_env_hilb_env = pmtm(env_corr_hilb_env, NW, NFFT, fs2);
P_env_hilb_env = pwelch(env_corr_hilb_env,WINDOW,OVERLAP,NFFT, fs2);
P_env_hilb_env = 10*log10(P_env_hilb_env);

%[s_2,c,~,~,~] = cmtm(input, phi,1/fs2,25); 

%% plot
close all 

simtime = length(psth_pos)/Fs;
tvect = 0:psthbinwidth:simtime-psthbinwidth;

tt= 0:1/Fs:(length(psth_pos)-1)/Fs;
t_env = 1/fs2:1/fs2:T;
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
plot(tt,px,t_env,hilb_env)
ylabel('Pressure (Pa)')
xlabel('Time (ms)')
title('AN Reponse to A4 - Violin | CF = 440 Hz')


figure;
ax1 = subplot(3,1,1);
%semilogx(freqPSTH, s_psd,'b');
%plotyy(freq_pmtm_env, sig_psd_env,freq_pmtm_tfs, sig_psd_tfs)
plot(freq_pmtm_sig, sig_psd_pmtm+100,'k',freq_pmtm_env, sig_psd_env,'b',freq_pmtm_tfs, sig_psd_tfs,'r')
set(gca, 'XScale', 'log')
title('Input PSD')
ylabel('PSD (dB/Hz)');
legend('Signal','Envelope','Fine Structure')
ylim([-200,0]);
%xlim([0,5000])

ax2 = subplot(3,1,2);
%semilogx(freqPSTH, s_psd,'b');
plot(freq_pmtm, s_psd_pmtm,'b')
set(gca, 'XScale', 'log')
title('Modulation PSD')
ylabel('PSD (dB/Hz)');
%xlim([0,5000])


ax3 = subplot(3,1,3);
%semilogx(freqPSTH, phi_psd,'r');
plot(freq_pmtm, phi_psd_pmtm,'r')
set(gca, 'XScale', 'log')
title('Carrier Freq PSD')
xlabel('Freq (Hz)');
ylabel('PSD (dB/Hz)');
%xlim([0,5000])

ax1.XLim = [1, 5000];
ax2.XLim = [1, 5000];
ax3.XLim = [1, 5000];
linkaxes([ax1,ax2,ax3],'x');

figure;
hold on 
semilogx(f_coherence,ENV_coherence, 'LineWidth',1.5);
semilogx(f_coherence,TFS_coherence, 'LineWidth',1.5);
hold off
set(gca, 'XScale', 'log')

xlabel('Frequency (Hz)')
ylabel('Coherence')

legend('ENV','TFS');

figure;
hold on 
semilogx(freq_pmtm,P_env_hilb_env,'LineWidth',1.5)
semilogx(freq_pmtm,P_tfs_hilb_tfs,'LineWidth',1.5);
hold off
set(gca, 'XScale', 'log')

xlabel('Frequency (Hz)')
ylabel('Cross-Spectral Density (dB/Hz)')
title('CSD of ENV/TFS Spectra with Hilbert ENV/TFS Spectra');
legend('ENV','TFS');
grid on
