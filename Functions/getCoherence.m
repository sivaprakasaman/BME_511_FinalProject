function [env_cohere, tfs_cohere, freq_cohere, s, phi, input_env, input_tfs] = getCoherence(bankedSig, input_fs, psth_pos, psth_neg, psth_fs, NFFT, CF);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    input_dur = length(bankedSig)/input_fs;
    psth_pos = psth_pos(:,1:ceil(input_dur*psth_fs));
    psth_neg = psth_neg(:,1:ceil(input_dur*psth_fs));

    input_env = abs(hilbert(bankedSig'));
    input_tfs = cos(angle(hilbert(bankedSig')));
    
    input_env = resample(input_env,psth_fs,input_fs);
    input_tfs = resample(input_tfs,psth_fs,input_fs);

    s = (psth_pos+psth_neg)./2;
    s = s';
    
    d = (psth_pos-psth_neg)./2;
    d = d';
    
    %questionable
    for i = 1:length(CF)
        [B,A] = butter(6,CF(i)/(psth_fs/2),'low');
        s(:,i) = filter(B,A,s(:,i));
    end
    
    phi = sqrt(2).*rms(d).*(d./abs(hilbert(d)));
    
    [env_cohere, freq_cohere] = mscohere(input_env,s,[],[],NFFT,psth_fs);
    [tfs_cohere, ~] = mscohere(input_tfs,d,[],[],NFFT,psth_fs);
    
    env_cohere = env_cohere.*~(freq_cohere>CF);
    
end

