
function [GNE] = GNE_measure(data, fs)

filt_order=100;
new_fs=10000;
x = 0.03*new_fs;
tstep=0.01*new_fs;
BW = 1000; % bandwidth
Fshift = 500; % shift fr   
data = resample(data, new_fs, fs); 

% cut-off1
Fc1=1:Fshift:(new_fs/2-BW-500); % not to cross Nq freq!
% cut-off2
Fc2=Fc1+BW;

for j=1:length(Fc1);
    d(j) =fdesign.bandpass('n,fc1,fc2',filt_order,Fc1(j),Fc2(j),new_fs);
    hd(j) = design(d(j));
end

steps=(length(data)-x)/tstep;

for i=1:steps+1   
    tseries = data(1+(i-1)*tstep:(i-1)*tstep+x);
    Dwindow = hann(length(tseries));
    segment_sig = tseries.*Dwindow;

    a = lpc(segment_sig,13);
    est_x = filter([0 -a(2:end)],1,segment_sig);    % Estimated signal
    e = segment_sig - est_x;
    LPE = xcorr(e,'coeff');   % LPES
    LPE=LPE(round(length(LPE)/2):end);

    for ii=1:length(hd)
        sigBW(:,ii)=filter(hd(ii),LPE);
        sig_TKEO(ii) = mean(TKEO(sigBW(:,ii)));
        sig_energy(ii) = mean(sigBW(:,ii)).^2; 
    end
    Hilb_tr = hilbert(sigBW);
    Hilb_env = abs(Hilb_tr);
    c = xcorr(Hilb_env);
    [cval,cidx] = max(c);
    GNEm(i) = max(cval);
    
    signal_BW_TKEO(i,:) = sig_TKEO;  
    signal_BW_energy(i,:) = sig_energy;
end

signal_BW_TKEO2 = mean(log(signal_BW_TKEO)); % used for getting the noise to signal ratio
signal_energy2 = mean(log(signal_BW_energy));

% Set outputs

GNE(1) = mean(GNEm);
GNE(2) = std(GNEm);

gnTKEO = mean(signal_BW_TKEO);
gnSEO = mean(signal_BW_energy);
GNE(3) = sum(gnTKEO(1:2))/sum(gnTKEO(end-3:end));
GNE(4) = sum(gnSEO(1:2))/sum(gnSEO(end-3:end));
GNE(5) = sum(signal_BW_TKEO2(end-3:end))/sum(signal_BW_TKEO2(1:2));
GNE(6) = sum(signal_energy2(end-3:end))/sum(signal_energy2(1:2));

end

