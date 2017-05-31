% tests count_ts2.m by running different cases to run through each of the 
% significant conditional statements in count_ts2.m

clear

% ts_tested = time shifts that the algorithm will test for
m = 2; ns = 10; 
for i = 1:ns % pick time shifts (treat each shift as a different spike type)
    ts(i,1:m) = 2*i-1; % time shifts
end

disp('------ CASE 1 (2 CORRECT, 1 FP, 1 CORRECT_TS, SAD = [-1 4]) -------------')
% goes through lines 76, 80, 90-91
[TN, FN, CORRECT, FP, WRONG, CORRECT_TS, SAD] = count_ts2([1 1], [1 11 4], ts, [2 3], 0)

disp('------ CASE 2 (2 CORRECT, 1 FP, 2 CORRECT_TS, SAD = [-1 1]) -------------')
% goes through lines 76, 80, 90-91
[TN, FN, CORRECT, FP, WRONG, CORRECT_TS, SAD] = count_ts2([1 1], [1 11 1], ts, [2 0], 0)

disp('------ CASE 3 (1 CORRECT, 1 FP, 1 WRONG, 1 CORRECT_TS, SAD = [-1]) ------')
% goes through lines 76, 80, 88-91
[TN, FN, CORRECT, FP, WRONG, CORRECT_TS, SAD] = count_ts2([1 1], [1 11 11], ts, [2 0], 0)

disp('------ CASE 4 (2 CORRECT, 1 FP, 2 CORRECT_TS, SAD = [-1 1]) -------------')
% goes through lines 76, 80, 88-91
[TN, FN, CORRECT, FP, WRONG, CORRECT_TS, SAD] = count_ts2([2 2], [1 11 11], ts, [2 0], 0)

disp('------ CASE 5 (1 CORRECT, 2 FP, SAD = [4) -------------------------------')
% goes through lines 76, 80, 88-91
[TN, FN, CORRECT, FP, WRONG, CORRECT_TS, SAD] = count_ts2(2, [1 11 17], ts, 9, 0)