%Code to inspect stimulus sound properties:

clear all, close all

addpath('Sound_Samples/Part A/')
addpath('Functions')

%for Part A:
%instruments = ["banjo","bassoon","cello","clarinet","flute","oboe","trumpet","saxophone","viola","violin"];
%instruments = ["violin","viola","cello"]; %String subset
instruments = ["banjo","clarinet","flute","trumpet","violin"];

pitch = 'A4';
cond = 'normal';
timewindow = [0,1.0];
nfft = 2048;

for i = 1:length(instruments)
    
    filename = strcat(instruments(i),'_',pitch,'_',cond,'.mp3')
    [DFTsig(i,:), DFTfreq_Hz, dataStruct, ~] = compute_dft(filename,timewindow(1),timewindow(2),nfft,'dB');
    fs = dataStruct.fs_Hz;
    sig = dataStruct.sig;
    
    hold on
    plot(DFTfreq_Hz,DFTsig(i,:));
    title(instruments(i));
    hold off
    sound(sig(timewindow(1)*fs+1:timewindow(2)*fs),fs);
    pause(1);
end

legend(instruments)
title("Frequency Response of A4 on Multiple Instruments");
