function y = synth(t, a, w, s, f1, f2, N, b)
% Makes signal vector (no noise) resulting from L events given firing times 
% (and other parameters), combines signals of overlapping spikes
%
% INPUTS: 
%         t       =  1 x L vector, firing times
%         a       =  1 x L vector, spike amplitudes
%         w       =  1 x L vector, spike width per firing event, 20*ceil(w)+1 samples
%         s       =  1 x L vector, spike shapes: 1 for f1, 0 for f2
%         f1, f2  =  functions @(t, w) describing spike shape
%         N       =  total number of sample times
%         b (opt) =  half the number of samples for spike width of 1
%
% OUTPUTS: 
%         y       =  1 x N vector, signal strength at each sample time

if nargin < 8, b = 10.0; end % default value of 10

y = zeros(1,N); % default signal strength is 0
for j = 1:length(t) % iterate through number of events   
    
    x = [floor(t(j)-b*w(j)):ceil(t(j)+b*w(j))]; % indices in output array y
    x(x<1) = [];  x(x>N) = []; % do not go beyond the sample time interval
    
    if s(j) == 1 % create spikes with shape f1
        y(x) = y(x) + a(j)*f1((x-t(j))/b,w(j)); 
    else % create spikes with shape f2
        y(x) = y(x) + a(j)*f2((x-t(j))/b,w(j)); 
end

end