function [TN, FN, CORRECT, FP, WRONG, CORRECT_TS, SAD] = count_ts2(spike_act, spike_det, ts_tested, ts_act, norm_frac)
% Calculates error fractions for spike type detection with time shifts
% (multiple spike types)
% 
% INPUTS: (a vector may also be inputed)
%       spike_act        = cell array of actual spike shapes
%       spike_det        = cell array of detected spike shapes and time
%                          shifts (counts each time shift as separate type), indices of ts_tested
%       ts_tested        = ns x m matrix of time shifts that the algorithm tests for
%       ts_act           = cell array of time shifts of the actual spikes
%       norm_frac (opt)  = 0 to keep raw count, 1 (default) to normalize
%                          error fractions
%
% OUTPUTS:
%        TN         = true negative, no spike present or detected
%           ^ normalized by dividing by (TN + FN) 
%
%        FN         = false negative, don't detect spike when one is present
%        CORRECT    = identify correct spike type
%        WRONG      = identify wrong spike type
%           ^ normalized by dividing by (FN + CORRECT + WRONG)
%
%        FP         = false positive, detect spike that isn't there
%           ^ normalized by dividing by (CORRECT + FP + WRONG)
% 
%        CORRECT_TS = identify the correct time shift (within 3 sample times from actual)
%           ^ normalized by dividing by total number of time shifts = total
%             number of actual spikes
%
%        SAD        = standard absolute deviation, best guess time shift - actual time shift 

if ~iscell(spike_act), spike_act = {spike_act}; end % if a vector is inputed, change to cell array
if ~iscell(spike_det), spike_det = {spike_det}; end
if ~iscell(ts_act), ts_act = {ts_act}; end
if nargin < 5, norm_frac = 1; end % default, normalize error fractions

[ns, m] = size(ts_tested); % number of time shifts, number of spike types
R = length(spike_act); % number of clips (runs)
TN = 0; FN = 0; CORRECT = 0; FP = 0; WRONG = 0; CORRECT_TS = 0; WRONG_TS = 0; % initialize spike counting
SAD = [];
num_act = 0; % start with no actual spikes
k = 1; % counter for # of actual time shifts
for i = 1:R % for each clip (run)
    act_sort = sort(spike_act{i}); % sorted list of actual spikes
    A = length(act_sort); % number of actual spikes
    num_act = num_act + A; % total # of actual spikes = total # of actual time shifts
    
    det_sort = sort(spike_det{i}); % sorted list of detected spikes 
    ts_det = ts_tested(det_sort); % detected time shifts
    det_sort = ceil(det_sort/ns); % convert to list of spike types only (don't include shifts as diff. types)
    D = length(det_sort); % number of detected spikes
    
    % counting error types for spike types
    if A == 0 && D == 0 % no spike present or detected = TN
        TN = TN + 1;
    elseif D == 0 && A > 0 % no spikes detected when spikes are present = FN
        FN = FN + A;
    elseif A == 0 && D > 0 % detect spikes when none are present = FP
        FP = FP + D;
    else % A,D are both nonzero
        a = 1; d = 1; % counters along act_sort, det_sort
        det_skip = 0; act_skip = 0; % count how many in each list is skipped over (i.e. no match is found)
        s = 0; % stop condition
        while ~s
            if a > A 
                if d <= D, det_skip = det_skip + D-d+1; end % skipped detected spikes
                s = 1; % stop checking, not enough spikes to compare
            elseif d > D 
                if a <= A, act_skip = act_skip + A-a+1; end % skipped actual spikes
                s = 1; % stop checking, not enough spikes to compare
            elseif det_sort(d) < act_sort(a)
                d = d + 1; % advance one in det_sort list
                det_skip = det_skip + 1; % no match found for detected spike type
            elseif det_sort(d) > act_sort(a)
                a = a + 1; % advance one in act_sort list
                act_skip = act_skip + 1; % no match found for actual spike type
            else % det_sort(d) == act_sort(a)
                CORRECT = CORRECT + 1; % match found
                [mm,nn] = min(abs(ts_det(det_sort == det_sort(d)) - ts_act{i}(a)));
                index = find(det_sort == det_sort(d)); % indices of detected time shifts that correspond with detected spike type
                if mm <= 3
                    CORRECT_TS = CORRECT_TS + 1;
                end
                SAD(k) = ts_det(index(nn)) - ts_act{i}(a); k = k + 1; % deviation from actual time shift
%                 if abs(SAD(k-1)) > 8
%                     disp('spike type'); disp(det_sort(d));
%                     disp('actual time shift'); disp(ts_act{i}(a)); disp(ts_act{i}); disp(det_sort)
%                     disp('detected time shift'); disp(ts_det(index(nn))); disp(ts_det)
%                 end
                ts_det(index(nn)) = NaN; % don't count this time shift again
                d = d + 1; a = a + 1; % advance one in both lists
            end
        end % end while loop
        w = min(det_skip,act_skip); % wrong identification will appear skipped on both lists
        WRONG = WRONG + w;
        if det_skip > act_skip % detect more than present
            FP = FP + (det_skip - act_skip);
        elseif det_skip < act_skip % detect fewer than present
            FN = FN + (act_skip - det_skip);
        end % end if statement
    end % end A,D conditions
end % end for loop over runs

if norm_frac == 0, % report raw count (NOT normalized)
    total_det = 1; total_emp = 1; total_act = 1; total_ts = 1;
else % normalize error fractions (default)
    total_emp = TN + FN; % total number of empty clips
    total_act = FN + CORRECT + WRONG; % total number of actual spikes
    total_det = CORRECT + FP + WRONG; % total number of detected spikes
    total_ts = num_act; % total number of time shifts = total number of actual spikes
end

% normalization of error fractions
if total_emp ~= 0 % don't divide by zero
    TN = TN/total_emp; % frac not detected that actually are not spikes
end
if total_det ~= 0 % don't divide by zero
    FP = FP/total_det; % frac detected that actually are not spikes
end
if total_act ~= 0 % don't divide by zero
    CORRECT = CORRECT/total_act; % frac actual spikes found + ID'ed correctly
    FN = FN/total_act; % frac actual spikes missed
    WRONG = WRONG/total_act; % frac actual spikes that are ID'ed incorrectly
end
if total_ts ~= 0 % don't divide by zero
    CORRECT_TS = CORRECT_TS/total_ts; % frac of correctly ID'ed time shifts
end
end