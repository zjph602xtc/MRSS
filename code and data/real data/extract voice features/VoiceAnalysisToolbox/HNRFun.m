
function [HNR, NHR] = HNRFun(data,fs)

f0max=500; %Hz -- max value, possibly adjust for other applications
f0min=50; %Hz
tstep=0.01*fs;
x=0.08*fs;
steps=(length(data)-x)/tstep;

for i=1:steps
    
    tseries = data(i*tstep:i*tstep+x);
    tseries = tseries-mean(tseries);
    Dwindow = hann(length(tseries));
    segment_sig = tseries.*Dwindow;
    
    %% HNR computation process
    ACF = xcorr(segment_sig,'coeff');
    ACF2 = ACF(length(segment_sig):end);
    aa=fft(segment_sig);
    aa=ifft(abs(aa).^2);
    ACF_Dwindow = xcorr(Dwindow,'coeff');
    ACF_Dwindow2 = ACF_Dwindow(length(Dwindow):end);
    bb=fft(Dwindow);
    bb=ifft(abs(bb).^2);
    ACF_signal = ACF2./ACF_Dwindow2;
    ACF_signal = ACF_signal(1:round(length(ACF_signal)/3));
    rho=aa./bb;
    rho=rho(1:round(length(rho)/2));
    rho=rho/max(rho);
    [rx_value,rx_index] = sort(ACF_signal,'descend');
    [d1 d2] = sort(rho, 'descend');
    low_lim=ceil(fs/f0max);  % round towards positive sample number
    up_lim=floor(fs/f0min);  % round towards negative sample number
    k=2;
    while ((rx_index(k)<low_lim) || rx_index(k)>up_lim)
        k=k+1;
    end
    
    m=2;
    while ((d2(m)<low_lim) || d2(m)>up_lim)
        m=m+1;
    end
    ll(i)=d2(m);
    mm=d2(m); 
    HNR_dB_Praat(i) = 10*log10(rho(mm)/(1-rho(mm)));
    NHR_Praat(i) = (1-rho(mm))/rho(mm);

end

%% Summarize data
HNR(1)=mean(HNR_dB_Praat);
HNR(2)=std(HNR_dB_Praat);

NHR(1)=mean(NHR_Praat);
NHR(2)=std(NHR_Praat);

end
