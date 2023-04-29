
function [GQ] = glottis_quotient(VF_close, VF_open, fs, f0min, f0max, flag)
%% Calculate the glottis quotients

cycle_open=abs(VF_open(2:end)-VF_close(1:end-1));
cycle_closed=abs(VF_open(1:end-1)-VF_close(1:end-1));

% remove erroneous cycles
if flag
    low_lim=fs/f0max; % lower limit
    up_lim=fs/f0min;  % upper limit
    N=length(cycle_open);
    for i=1:N-1
        if((cycle_open(i) > up_lim) || (cycle_open(i) < low_lim))
            cycle_open(i)=NaN;
        end
        if((cycle_closed(i) > up_lim) || (cycle_closed(i) < low_lim))
            cycle_closed(i)=NaN;
        end
    end
end

%statistics in time
prc1=prctile(cycle_open,[5 95]);
cycle_open_range_5_95_perc=prc1(2)-prc1(1);
prc2=prctile(cycle_closed,[5 95]);
cycle_closed_range_5_95_perc=prc2(2)-prc2(1);

GQ(1) = (cycle_open_range_5_95_perc/(cycle_open_range_5_95_perc+cycle_closed_range_5_95_perc));
GQ(2) = (nanstd(cycle_open));
GQ(3) = (nanstd(cycle_closed));

end
