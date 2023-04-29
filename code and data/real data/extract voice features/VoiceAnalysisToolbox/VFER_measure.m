
function [VFER] = VFER_measure(data, fs)

filt_order=100;
BW = 500; % bandwidth 
Fmax = (fs/2-BW-300); % Max frequency to check 
Fshift = 500; % shift 

% get VF action
[VF_close,VF_open] = dypsa(data,fs);

% cut-off1
Fc1=1:Fshift:Fmax;
% cut-off2
Fc2=Fc1+BW;

for j=1:length(Fc1);
    d(j) = fdesign.bandpass('n,fc1,fc2',filt_order,Fc1(j),Fc2(j),fs);
    hd(j) = design(d(j));
end

for i=1:length(VF_close)-1
    tseries = data(VF_close(i):VF_close(i+1));
    Dwindow = hann(length(tseries)); %Use Hanning window
    segment_sig = tseries.*Dwindow;
    
    if (length(tseries)>50)
        for ii=1:length(hd)
            thanasis = filter(hd(ii),segment_sig);
            sigBW(:,ii) = thanasis(1:50);
            sig_TKEO(ii) = mean(TKEO(sigBW(:,ii)));
            sig_SEO(ii) = mean(sigBW(:,ii)).^2; 
        end
        Hilb_tr = hilbert(sigBW);
        Hilb_env = abs(Hilb_tr);
        c = xcorr(Hilb_env);
        [cval,cidx] = max(c);
        NEm(i) = max(cval);

        signal_BW_TKEO(i,:) = sig_TKEO;
        signal_BW_SEO(i,:) = sig_SEO;        
    end
end

signal_BW_TKEO2 = mean(log(signal_BW_TKEO)); % used for getting the noise to signal ratio

% Set outputs

VFER(1) = mean(NEm);
VFER(2) = std(NEm);
VFER(3) = -sum(NEm.*log_bb(NEm));
VFTKEO = mean(signal_BW_TKEO);
VFSEO = mean(signal_BW_SEO);
VFlog_SEO = mean(log(signal_BW_SEO));

% Get 'signal to noise' ratios
VFER(4) = sum(VFTKEO(1:5))/sum(VFTKEO(6:10));
VFER(5) = sum(VFSEO(1:5))/sum(VFSEO(6:10));
VFER(6) = sum(signal_BW_TKEO2(6:10))/sum(signal_BW_TKEO2(1:5));
VFER(7) = sum(VFlog_SEO(6:10))/sum(VFlog_SEO(1:5));

end
