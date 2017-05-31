% tests count.m by running different cases to run through each of the 
% significant conditional statements in count.m

clear

disp('------ CASE 1 (1 TN) ------') % goes through line 37
[TN, FN, CORRECT, FP, WRONG] = count([],[],0)

disp('------ CASE 2 (2 FN) ------') % goes through line 39
[TN, FN, CORRECT, FP, WRONG] = count([1 4],[],0)

disp('------ CASE 3 (1 FP) ------') % goes through line 41
[TN, FN, CORRECT, FP, WRONG] = count([],4,0)

disp('------ CASE 4 (2 C, 1 FP, 1 W) ------') % goes through lines 54, 60, 67
[TN, FN, CORRECT, FP, WRONG] = count([3 3 3],[1 2 3 3],0)

disp('------ CASE 5 (1 C, 2 FP, 1 W) ------') % goes through lines 48-9, 54, 57, 60, 67-8
[TN, FN, CORRECT, FP, WRONG] = count([3 3],[4 1 2 3],0)

disp('------ CASE 6 (2 C, 2 FN) ------') % goes through lines 51-2, 57, 60, 70
[TN, FN, CORRECT, FP, WRONG] = count([1 3 3 4],[3 3],0)
