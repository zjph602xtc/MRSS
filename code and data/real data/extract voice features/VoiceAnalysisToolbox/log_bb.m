
function pout = log_bb(pin, method)
% Function that computes the algorithm depending on the user specified
% base; if the input probability is zero it returns zero.

if nargin<2
    method = 'Nats';
end

switch (method)
    case 'Hartson' % using log10 for the entropy computation
        log_b=@log10;
        
    case 'Nats' % using ln (natural log) for the entropy computation 
        log_b=@log;
       
    otherwise % method -> 'Bits' using log2 for the entropy computation 
        log_b=@log2;
end

if pin==0
    pout=0;
else
    pout=log_b(pin);
end

end