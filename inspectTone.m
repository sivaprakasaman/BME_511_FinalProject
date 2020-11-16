%Code to inspect stimulus sound properties:

clear all, close all

addpath('Sound_Samples/Part A/')
addpath('Sound_Samples/Part B/Tambourine')

addpath('Functions')

%for Part A:
%instruments = ["banjo","bassoon","cello","clarinet","flute","oboe","trumpet","saxophone","viola","violin"];
%instruments = ["violin","viola","cello"]; %String subset
instruments = ["banjo","clarinet","flute","trombone","violin"];
%instruments = ["tambourine"];
pitch = 'A4';
cond = 'normal';
%cond = "long_cresc-decresc_shaken";
timewindow = [0,1.0];
nfft = 2048;

for i = 1:length(instruments)
    
    filename = strcat(instruments(i),'_',pitch,'_',cond,'.mp3')
    [DFTsig(i,:), DFTfreq_Hz, dataStruct, ~] = compute_dft(filename,timewindow(1),timewindow(2),nfft,'mag');
    fs = dataStruct.fs_Hz;
    sig = dataStruct.sig;
    sig = sig(timewindow(1)*fs+1:timewindow(2)*fs);
    hold on
    plot(DFTfreq_Hz,DFTsig(i,:));
    hold off
    sound(sig,fs);
%     getSpect(sig,40,fs,65,'mag',strcat(instruments(i),' ',' Crescendo_Hand'))
%     cd Figures
%     saveas(gcf,strcat('spectrogram_',instruments(i),'_',cond),'epsc')
%     cd ../
%     close all;
    pause(1);
    
end

legend(instruments)
title("Frequency Response of A4 on Multiple Instruments");
xlim([0 max(DFTfreq_Hz)]);
xlabel("Frequency (Hz)")
ylabel("Amplitude")
