
function [IMF] = IMF_measure(data)

%% Use classical EMD

IMF_dec = emd(data);
IMF_dec=IMF_dec';
IMF_dec2=log_bb(IMF_dec); %Log transformation
[N,M]=size(IMF_dec);

for i=1:M
    IMF_decEnergy(i) = abs(mean((IMF_dec(:,i)).^2));
    IMF_decTKEO(i) = abs(mean(TKEO(IMF_dec(:,i))));
    IMF_decEntropia(i) = abs(mean(-sum(IMF_dec(:,i).*log_bb(IMF_dec(:,i)))));
    IMF_decEnergy2(i) = abs(mean((IMF_dec2(:,i)).^2));
    IMF_decTKEO2(i) = abs(mean(TKEO(IMF_dec2(:,i))));
    IMF_decEntropia2(i) = abs(mean(-sum(IMF_dec2(:,i).*log_bb(IMF_dec2(:,i)))));
end
    
% Get 'signal to noise' ratio measures
IMF(1) = sum(IMF_decEnergy(4:end))/sum(IMF_decEnergy(1:3));
IMF(2) = sum(IMF_decTKEO(4:end))/sum(IMF_decTKEO(1:3));
IMF(3) = sum(IMF_decEntropia(4:end))/sum(IMF_decEntropia(1:3));
IMF(4) = abs(sum(IMF_decEnergy2(1:2))/sum(IMF_decEnergy2(4:end)));
IMF(5) = abs(sum(IMF_decTKEO2(1:2))/sum(IMF_decTKEO2(3:end)));
IMF(6) = sum(IMF_decEntropia2(1:2))/sum(IMF_decEntropia2(3:end));

end
