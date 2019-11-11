function stats = rm_anova2_matrix(Y,S,F1,F2)
%
% function stats = rm_anova2(Y,S,F1,F2,FACTNAMES)
%
% Two-factor, within-subject repeated measures ANOVA.
% For designs with two within-subject factors.
%
% Parameters:
%    Y          dependent variable (numeric) in a column vector
%    S          grouping variable for SUBJECT
%    F1         grouping variable for factor #1
%    F2         grouping variable for factor #2
%    FACTNAMES  a cell array w/ two char arrays: {'factor1', 'factor2'}
%
%    Y should be a 1-d column vector with all of your data (numeric).
%    The grouping variables should also be 1-d numeric, each with same
%    length as Y. Each entry in each of the grouping vectors indicates the
%    level # (or subject #) of the corresponding entry in Y.
%
% Returns:
%    stats is a cell array with the usual ANOVA table:
%      Source / ss / df / ms / F / p
%
% Notes:
%    Program does not do any input validation, so it is up to you to make
%    sure that you have passed in the parameters in the correct form:
%
%       Y, S, F1, and F2 must be numeric vectors all of the same length.
%
%       There must be at least one value in Y for each possible combination
%       of S, F1, and F2 (i.e. there must be at least one measurement per
%       subject per condition).
%
%       If there is more than one measurement per subject X condition, then
%       the program will take the mean of those measurements.
%
% Aaron Schurger (2005.02.04)
%   Derived from Keppel & Wickens (2004) "Design and Analysis" ch. 18
%

%
% Revision history...
%
% 11 December 2009 (Aaron Schurger)
% 
% Fixed error under "bracket terms"
% was: expY = sum(Y.^2);
% now: expY = sum(sum(sum(MEANS.^2)));
%
% 05 April 2019 (Roger Strong)
% Can input multiple columns of Y at once (do separate ANOVA on each
% column, useful for simulations)
% Have a few slow steps (creating AS and BS with for loops)

F1_lvls = unique(F1);
F2_lvls = unique(F2);
Subjs = unique(S);

a = length(F1_lvls); % # of levels in factor 1
b = length(F2_lvls); % # of levels in factor 2
n = length(Subjs); % # of subjects

AB = zeros(a,b,size(Y,2));

for i = 1:a
    for j = 1:b
        AB(i, j, :) = sum(Y(F1 == i & F2 == j, :));
    end
end
    

%original
AS = zeros(a, n, size(Y,2));
for i = 1:a
    for sub = 1:n
        AS(i, sub, :) = sum(Y(S == sub & F1 == i, :));
    end
end
BS = zeros(b, n, size(Y,2));
for j = 1:b
    for sub = 1:n
        BS(j, sub, :) = sum(Y(S == sub & F2 == j, :));
    end 
end

%%%tried without subject (n) for loop, even slower though
% %new
% Y2 = reshape(Y,numel(Y),1);
% S2 = repmat(S, size(Y,2), 1);
% F12 = repmat(F1, size(Y,2), 1);
% F22 = repmat(F2, size(Y,2), 1);
% 
% AS2 = zeros(a, n, size(Y,2));
% for i = 1:a
%     tmp = Y2(F12 == i);
%     AS2(i, :, :) = sum(permute(reshape(tmp, n, a, size(Y,2)), [2 1 3]));
% end
% BS2 = zeros(b, n, size(Y,2));
% for i = 1:b
%     tmp = Y2(F22 == i);
%     BS2(i, :, :) = sum(permute(reshape(tmp, n, b, size(Y,2)), [2 1 3]));
% end


A = squeeze(sum(AB,2)); % sum across columns, so result is ax1 column vector
B = permute(squeeze(sum(AB,1)), [2 1]); % sum across rows, so result is 1xb row vector
S = permute(squeeze(sum(AS,1)), [2 1]); % sum across columns, so result is 1xs row vector
T = sum(A); % could sum either A or B or S, choice is arbitrary

% degrees of freedom
dfA = a-1;
dfB = b-1;
dfAB = (a-1)*(b-1);
dfS = n-1;
dfAS = (a-1)*(n-1);
dfBS = (b-1)*(n-1);
dfABS = (a-1)*(b-1)*(n-1);

% bracket terms (expected value)
expA = (sum(A.^2)./(b*n))';
expB = sum(B.^2, 2)./(a*n);
expAB = squeeze(sum(sum(AB.^2),2)./n);


expS = sum(S.^2, 2)./(a*b);
expAS = squeeze(sum(sum(AS.^2), 2)./b);
expBS = squeeze(sum(sum(BS.^2), 2)./a);
expY = sum(Y.^2)';
expT = (T.^2 ./ (a*b*n))';

% sums of squares
ssA = expA - expT;
ssB = expB - expT;
ssAB = expAB - expA - expB + expT;
ssS = expS - expT;
ssAS = expAS - expA - expS + expT;
ssBS = expBS - expB - expS + expT;
ssABS = expY - expAB - expAS - expBS + expA + expB + expS - expT;
ssTot = expY - expT;

% mean squares
msA = ssA ./ dfA;
msB = ssB ./ dfB;
msAB = ssAB ./ dfAB;
msS = ssS ./ dfS;
msAS = ssAS ./ dfAS;
msBS = ssBS ./ dfBS;
msABS = ssABS ./ dfABS;

% f statistic
fA = msA ./ msAS;
fB = msB ./ msBS;
fAB = msAB ./ msABS;

% p values
stats.pA = 1-fcdf(fA,dfA,dfAS);
stats.pB = 1-fcdf(fB,dfB,dfBS);
stats.pAB = 1-fcdf(fAB,dfAB,dfABS);

% return values
% stats = {'Source','SS','df','MS','F','p';...
%          FACTNAMES{1}, ssA, dfA, msA, fA, pA;...
%          FACTNAMES{2}, ssB, dfB, msB, fB, pB;...
%          [FACTNAMES{1} ' x ' FACTNAMES{2}], ssAB, dfAB, msAB, fAB, pAB;...
%          [FACTNAMES{1} ' x Subj'], ssAS, dfAS, msAS, [], [];...
%          [FACTNAMES{2} ' x Subj'], ssBS, dfBS, msBS, [], [];...
%          [FACTNAMES{1} ' x ' FACTNAMES{2} ' x Subj'], ssABS, dfABS, msABS, [], []};
 
 return