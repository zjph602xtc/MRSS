
function PQ = perq1(time_series, K)

%% Calculate the PQ using the classical PQ formula

N = length(time_series);
mean_tseries = mean(time_series);
K1=round(K/2);
K2=K-K1;
p = 5;
sum1=0;

for i = K1:N-K2
    sum1 = sum1+mean(abs([time_series(i-K2:i+K2)]-time_series(i)));
end
        
PQ.classical_Schoentgen = (sum1/(N-K+1))/(mean_tseries);

sum2=0;
for i = K1:N-K2
    sum2 = sum2+mean(abs([time_series(i-K2:i+K2)]))-time_series(i);
end
        
PQ.classical_Baken = (sum2/(N-K+1))/(mean_tseries);

% perturbation quotient of the residue
time_series=time_series(:);
sum3=0;
% calculate the AR coefficients (I use the Yule-Walker equations)
new_tseries=(time_series-mean_tseries)';
a = aryule(time_series-mean_tseries,p);

for i = 1+p:N
    sum3 = sum3+abs(sum(a.*(new_tseries(i:-1:i-p))));
end

PQ.generalised_Schoentgen = (sum3/(N-p))/(mean_tseries);
    
end
