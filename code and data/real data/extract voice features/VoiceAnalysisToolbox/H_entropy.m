
function H = H_entropy(f)
% Calculate entropy of a discrete distribution
% Usage: H = entropy(f)
%  f - input distribution as a vector
%  H - entropy
N = length(f);
H = 0;
for j = 1:N
   H = H - f(j) * log_bb(f(j));
end

end % end of H_entropy function
