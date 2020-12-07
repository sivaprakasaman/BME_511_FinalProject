function [] = getSpect(sig,BW,Fs,DR,spectType,sigName)
%Modified version of spectrogram_BW_DR.m from Heinz
%Added an FFT for better visualization of peaks 

%TODO: CLEAN UP!!

% Input Parameters
%   sig: waveform to be analyzed
%   BW : bandwidth of spectrogram filters
%   Fs: sampling rate of signal
%   DR: limits the dynamic range of Sgram to improve visualization of amplitude differences
%
% Output Parameters
%   Sgram: matrix representation of spectrogram
%          - each column is the Short-time Fourier Transform (STFT) at a given point in time (energy vs frequency) 
%          - each row represents the energy at a given frequency as a function of time
%   SG_Freq_Hz  : Frequency vector corresponding to rows in Sgram
%   SG_Time_sec : Time vector corresponding to columsn in Sgram



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transform the BW_Hz parameter into the parameters MATLAB's "spectrogram" 
% function needs, i.e., number of points.  And, specify default values for
% other parameters, which can be changed if needed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uses Hamming window as default.  This calculation is based on the
% fact that a Hamming window has a BW 1/2 as large as a rectangular window
% (which has a BW = 2/T, where T is the duration of the window)
Nwindow = round(4*Fs/BW);   % Number of points in Hamming window applied to signal prior to STFT calculation
% Set the default overlap of temporal windows to be 90%, i.e., each time
% step in Sgram is 10% if the temporal window length
PercentOverlap = 90;
OverlapFactor = (100-PercentOverlap)/100 + 1;
Noverlap = round(Nwindow/OverlapFactor);  % Number of points of overlap
% Set the number of points in each FFT, this is the default MATLAB uses
Nfft = max([256 2^(floor(log2(Nwindow))+1)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call MATLAB's spectrogram command with their parameters
% Call compute_dft to get an FFT as well 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Sgram,SG_Freq_Hz,SG_Time_sec] = spectrogram(sig,Nwindow,Noverlap,Nfft,Fs); 

nfft_spect = max([2048, 2^(floor(log2(length(sig))+1))]);
[DFTsig, DFTfreq_Hz, dataStruct, ~] = compute_dft(sig,Fs,[],[],nfft_spect,spectType);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the spectrogram and signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AnnotationFontSize = 11;
LabelFontSize = 11;

% Create the time_vector for the signal
Npts=length(sig);
time_signal_sec = (0:Npts-1)/Fs;


%% LIMIT DYNAMIC RANGE to improve visualization
SgramFactor = 10^-(DR/20);  % finds factor corresponding to DynamicRange_dB down
SgramFloor = SgramFactor*max(max(abs(Sgram)));  % finds value of Sgram corresponding to DynamicRange_dB down from peak
Sgram_dB = 20*log10(abs(Sgram) + SgramFloor);  % Creates a "noise floor" at DynamicRange_dB down

if ishandle(1),     close(1), end
figure(1); clf
%% This helps with the positioning and sizing of the figure 
set(gcf,'PaperPositionMode','auto')
set(gcf,'units','inches')
fig_xcorner_inches = 20;
fig_ycorner_inches = 10;
fig_xlength_inches = 7;
fig_ylength_inches = 8;
set(gcf,'Position',[fig_xcorner_inches fig_ycorner_inches fig_xlength_inches fig_ylength_inches]);

h1 = subplot(4,1,[1 2]);
surf(SG_Time_sec,SG_Freq_Hz,Sgram_dB,'EdgeColor','none');  % These commands are from MATLAB's help on "spectrogram"
axis xy; axis tight; colormap(parula); view(0,90);
% title(sprintf('Spectrogram [BW = %.f Hz; Dynamic Range = %.f dB]',BW,DR),'FontSize',LabelFontSize)
title(['Spectrogram:',sigName]);
xlabel('Time (sec)','FontSize',LabelFontSize);
ylabel('Frequency (Hz)','FontSize',LabelFontSize);
ylim([0 5000])
xlim([0 max(time_signal_sec)])
set(h1,'FontSize',AnnotationFontSize)  % set font size
%% Manipulate the position of FIG1 to use space efficiently
set(h1,'Units','norm')  % use normalized units
h1pos=get(h1,'pos');
% Yexpansion = 1.6;    % factor to stretch FIG1 vertically
% h1pos(2) = h1pos(2) - (Yexpansion-1) * h1pos(4);  % shift FIG1 down
% h1pos(4) = Yexpansion * h1pos(4);  % expand FIG1 vertically
% set(h1,'pos',h1pos)  % set the new FIG1 position

% figure(2);clf
h2 = subplot(413);
plot(time_signal_sec,sig)
ylabel('Amplitude','FontSize',LabelFontSize)
xlabel('Time (sec)','FontSize',LabelFontSize)
xlim([0 max(time_signal_sec)])
title('Signal waveform','FontSize',LabelFontSize)
set(h2,'FontSize',AnnotationFontSize)
linkaxes([h1,h2],'x')
format lonG
%plot spectrum:
h3 = subplot(414);
plot(DFTfreq_Hz,DFTsig)
ylabel('Amplitude','FontSize',LabelFontSize)
xlabel('Frequency (Hz)','FontSize',LabelFontSize)
xlim([0 5000])
xticks(440:440:max(DFTfreq_Hz));
title('FFT','FontSize',LabelFontSize)
set(h3,'FontSize',AnnotationFontSize)

ax = ancestor(h3, 'axes');
ax.XAxis.Exponent = 0;
%% Manipulate the position of FIG2 to match width of FIG1
set(h2,'Units','norm')  % use normalized units
h2pos=get(h2,'pos');
h2pos(2) = h2pos(2)-.03;
set(h3,'Units','norm')  % use normalized units
h3pos=get(h3,'pos');
h3pos(2) = h3pos(2)-.045;
%  h2pos(1) = h1pos(1);  % set FIG1 xcorner equal to FIG1 
% % h2pos(3) = h1pos(3);  % set FIG1 xlength equal to FIG1 
% % % h2pos(2)=.08;
set(h2,'pos',h2pos)  % set the new FIG2 position
set(h3,'pos',h3pos)  % set the new FIG3 position

% Listen to the sound
%sound(sig,Fs)

end

