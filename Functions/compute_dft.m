function [DFTsig, DFTfreq_Hz, dataStruct, paramsStruct] = compute_dft(varargin)

% Check if input data is in vector format or filename format 
    
    isFile = isstring(varargin{1})||ischar(varargin{1});
    
% If input data is in filename format, then 
%   - Read audiofile 
%   - Read or initialize tStart_sec, tEnd_sec, nfft, FFTtype
% Elseif input data is in vector format, then  
%   - Read or initialize data, fs_Hz, tStart, tEnd, nfft, FFTtype

if isFile
    
    switch nargin
        case 1
            [sig,fs_Hz] = audioread(varargin{1});
            tStart_sec = [];
            tEnd_sec = [];
            nfft = [];
            FFTtype = [];
            
        case 5
            [sig,fs_Hz] = audioread(varargin{1});
            tStart_sec = varargin{2};
            tEnd_sec = varargin{3};
            nfft = varargin{4};
            FFTtype = varargin{5};
    end
    
elseif isvector(varargin{1}) && nargin == 6        
  
        sig = varargin{1};
        fs_Hz = varargin{2};
        tStart_sec = varargin{3};
        tEnd_sec = varargin{4};
        nfft = varargin{5};
        FFTtype = varargin{6};
else
    disp('invalid inputs, try again');
    return
end 

% 
%     if isFile
%        
%         [sig,fs_Hz] = audioread(varargin{1});
%         tStart_sec = varargin{2};
%         tEnd_sec = varargin{3};
%         nfft = varargin{4};
%         FFTtype = varargin{5};
%     else
%         sig = varargin{1}; 
%         fs_Hz = varargin{2};
%         tStart_sec = varargin{3};
%         tEnd_sec = varargin{4};
%         nfft = varargin{5};
%         FFTtype = varargin{6};
%     end

    
    %Init Defaults (if not provided)
    t_total = length(sig)/fs_Hz;
    
    if(isempty(tStart_sec))
        tStart_sec = 0;
    end
    
    
    if(isempty(tEnd_sec))
        tEnd_sec = t_total;
    end
    
    
    if(isempty(FFTtype))
        FFTtype = 'mag';
    end
    
% Assign fields to output structure #3 (dataStruct) here
    dataStruct.sig = sig(:);
    dataStruct.fs_Hz = fs_Hz;
    
% Check if the values corresponding to tStart and tEnd are valid
    if ((tStart_sec >= tEnd_sec)||(tStart_sec<0)||(tEnd_sec > t_total))
       
        DFTsig = NaN;
        DFTfreq_Hz = NaN; 
        paramsStruct = NaN;
        
        return
    end
% Initialize time vector (first sample should start at time 0)
    %t = 0:1:(t_total)*fs_Hz;
% Find indices in t_inSig that are within tStart and tEnd. (Inlclude both end points)
    
    t_diff = tEnd_sec-tStart_sec;
    
    %Number of samples from start/end
    t_diff_vect = 1:1:t_diff*fs_Hz;
    
    %Time-shift samples by t_Start;
    t_outSig = t_diff_vect + tStart_sec*fs_Hz;

% Extract desired segment from signal 
    outSig = sig(t_outSig);
% Check if nfft is large enough to avoid aliasing, else set its value 

    if(isempty(nfft)||nfft<length(t_diff_vect))
        nfft = 2^nextpow2(length(t_diff_vect));
        warning('Use higher nfft to avoid aliasing!!');
    end

% Assign fields to output structure #4 (paramsStruct) here
    paramsStruct.tStart_sec = tStart_sec;
    paramsStruct.tEnd_sec = tEnd_sec;
    paramsStruct.nfft = nfft;
    paramsStruct.FFTtype = FFTtype;

% Get fft of the segment: call it DFTsig
    DFTsig = fft(outSig, nfft);

% Initialize frequency vector for DFTsig: call it DFTfreq_Hz
    %+
    DFTfreq_Hz = (fs_Hz/2)*linspace(0,1,ceil(nfft/2))';

% Keep only nonnegative frequencies for DFTsig and DFTfreq_Hz
    DFTsig = DFTsig(1:ceil(nfft/2));
% Transform DFTsig based on FFTtype
    
    if strcmp(FFTtype,'mag')
        
       DFTsig = abs(DFTsig);
        
    elseif strcmp(FFTtype,'dB')
        
       DFTsig = 20*log10(abs(DFTsig));
        
    end
    
 DFTsig = DFTsig(:);
 DFTfreq_Hz = DFTfreq_Hz(:);
end 