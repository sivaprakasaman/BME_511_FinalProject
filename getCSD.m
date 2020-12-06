function [env_csd, tfs_csd, env_psd, tfs_psd, input_env_psd, input_tfs_psd, freq] = getCSD(input_env, input_tfs, s, phi, NW, NFFT, psth_fs, estimate_psd)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    WINDOW = hamming(ceil(length(input_env)/3));
    OVERLAP = .75;

    if strcmp(estimate_psd,"pwelch")

        [input_env_psd, freq] =  pwelch(input_env,WINDOW,OVERLAP,NFFT,psth_fs);
        input_tfs_psd =  pwelch(input_tfs,WINDOW,OVERLAP,NFFT,psth_fs);

        env_psd =  pwelch(s,WINDOW,OVERLAP,NFFT, psth_fs);
        tfs_psd =  pwelch(phi,WINDOW,OVERLAP,NFFT, psth_fs);

    elseif strcmp(estimate_psd, "pmtm")
        [input_env_psd, freq] =  pmtm(input_env,NW, NFFT, psth_fs);
        input_tfs_psd =  pmtm(input_tfs,NW, NFFT, psth_fs);

        env_psd =  pmtm(s,NW, NFFT, psth_fs);
        tfs_psd =  pmtm(phi,NW, NFFT, psth_fs);
    end
    
    env_corr = zeros([2*size(input_env,1)-1,size(input_env,2)]);
    tfs_corr = zeros([2*size(input_env,1)-1,2*size(input_env,2)]);

    for i = 1:size(input_env,2)
        
        env_corr(:,i) = xcorr(input_env(:,i), s(:,i));
        tfs_corr(:,i) = xcorr(input_tfs(:,i), phi(:,i));
        
    end
    
    env_csd = pwelch(env_corr, WINDOW, OVERLAP, NFFT, psth_fs);
    tfs_csd = pwelch(tfs_corr, WINDOW, OVERLAP, NFFT, psth_fs);

    input_env_psd = 10*log10(input_env_psd);
    input_tfs_psd = 10*log10(input_tfs_psd);

    env_psd = 10*log10(env_psd);
    tfs_psd = 10*log10(tfs_psd);
    
    env_csd = 10*log10(env_csd);
    tfs_csd = 10*log10(tfs_csd);
    
end

