addpath('Sound_Samples/Part A/')
addpath('Chimera code and WAV files')
addpath('BEZ2018model')

%Shouldn't have to run this more than once right?

cd BEZ2018model
mexANmodel
cd ../

clear all; clc
close all

[test_sig test_fs]= audioread('violin_A4_normal.mp3');

%% Following the general structure of testANmodel_BEZ2018.m:

%TODO: Read more on these params
% model parameters
CF    = 440;   % CF in Hz; %gonna need to change this
spont = 70;   % spontaneous firing rate %SATYA CHANGED TO 70/s
tabs   = 0.6e-3; % Absolute refractory period
trel   = 0.6e-3; % Baseline mean relative refractory period
cohc  = 1.0;    % normal ohc function
cihc  = 1.0;    % normal ihc function
species = 2;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning) %read up on this tuning
noiseType = 0;  % 1 for variable fGn; 0 for fixed (frozen) fGn
implnt = 0;     % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

stimdb = 60; % stimulus intensity in dB SPL
F0 = CF;     % stimulus frequency in Hz
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 1;  % stimulus duration in seconds
rt = 2.5e-3; % rise/fall time in seconds
%ondelay = 10e-3;

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

Psth_pos = sum(reshape(psth_pos,psthbins,length(psth_pos)/psthbins)); %
Psth_neg = sum(reshape(psth_neg,psthbins,length(psth_neg)/psthbins)); %


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