addpath('Sound_Samples/')
addpath('Chimera code and WAV files')
addpath('BEZ2018model')

%Shouldn't have to run this more than once right?

cd BEZ2018model
mexANmodel
cd ../

clear all; clc

%% Following the general structure of testANmodel_BEZ2018.m:

%TODO: Read more on these params
% model parameters
CF    = 1e3;   % CF in Hz; %gonna need to change this
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
T  = 1.75;  % stimulus duration in seconds
rt = 2.5e-3; % rise/fall time in seconds
ondelay = 10e-3;

%Modify to work with alternating polarities

nrep = 25; % number of stimulus repetitions (e.g., 50); 
psthbinwidth = 1e-4; % binwidth in seconds;
psthbins = round(psthbinwidth*Fs);  % number of psth bins per psth bin

