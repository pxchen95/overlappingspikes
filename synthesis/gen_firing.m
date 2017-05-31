function t = gen_firing(N, L0, R)
% Generates firing times using a Poisson distribution
% 
% INPUTS: 
%    N          = total number of sample times
%    L0         = mean # of firing times in N samples
%    R          = refractory period, given in # of sample times
%
% OUTPUTS: 
%    t          = 1 x L vector of firing times, in ascending order

t(1) = rand*10; % randomly pick first firing time
L = 1; % counter, number of firing events
lambda = L0/(N); % mean firing rate

while t(end) < N
    L = L + 1; % update counter
    t_R = poissrnd(R + 1/lambda); % next firing event will occur a Poisson random time after refractory period
    t(L) = t(L-1) + t_R; % create next firing time
end

if t(end) > N, t(end) = []; end % remove out of range spikes

end

