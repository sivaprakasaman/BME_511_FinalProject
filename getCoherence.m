function [env_cohere, tfs_cohere, freq_cohere, s, phi, input_env, input_tfs] = getCoherence(bankedSig, input_fs, psth_pos, psth_neg, psth_fs, NFFT, CF);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
%     env_cohere = zeros(length(CF), NFFT/2+1);
%     tfs_cohere = zeros(length(CF), NFFT/2+1);
%     s = zeros(length(CF), dur*psth_fs);
%     phi = zeros(length(CF), dur*psth_fs);

    input_env = abs(hilbert(bankedSig'));
    input_tfs = cos(angle(hilbert(bankedSig')));
    
    input_env = resample(input_env,psth_fs,input_fs);
    input_env = resample(input_env,psth_fs,input_fs);

    s = (psth_pos+psth_neg)./2;
    d = (psth_pos-psth_neg)./2;
    phi = sqrt(2)*rms(d)*(d./abs(hilbert(d)));

    [env_cohere, freq_cohere] = mscohere(input_env',s',[],[],NFFT,psth_fs);
    [tfs_cohere, ~] = mscohere(input_tfs',phi',[],[],NFFT,psth_fs);
    

end

