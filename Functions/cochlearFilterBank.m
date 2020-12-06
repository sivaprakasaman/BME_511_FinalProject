function [filtered_sig] = cochlearFilterBank(input,Fs_Hz,CF,hrect)
%cochlearFilterbank
%   Returns a set of filtered vectors corresponding to a particular CF 
%   gammatone_c from Ning Ma: 
%    https://staffwww.dcs.shef.ac.uk/people/N.Ma/resources/gammatone/
%   [filtered_sig] = cochlearFilterBank(input,CF)
    
    if nargin<4
       hrect = 0;
    end
    
    filtered_sig = zeros(length(CF),length(input));
    envelopes = zeros(length(CF),length(input));

    for i = 1:length(CF)
        [filtered_sig(i,:), ~]= gammatone_c(input, Fs_Hz, CF(i), hrect);
    end
    

end

