function x = qpFeatureSelection(classes, features)

% Quadratic Programming Algorithms for small nubmer of varibale (M)
% When M is small number, QP is solved without approximation method.
% when M is large, Q should be approximated usgin its eigenvalue of decomposition
% 
% Rodriguez-lujan, I. Quadratic Programming Feature Selection. J. Mach. Learn. Res. 11, 1491?1516 (2010).
%
% Code provided by Tsuyoshi Mikkaichi <tsuyoshimikkaichi@gmail.com>
% and Alexander Hoffmann <ahoffmann@ucla.edu>
% adapted to function by Joao Xavier <xavierj@mskcc.org>
% April 10, 2019

n = length(unique(classes));
M = size(features, 2);

% correlation between feature and class
[r, pp]= corr(dummyvar(categorical(classes)),features);
%[r(2, :), pp]= corr(1-classes,features);

% correlations between features
[r2,p2]= corr(features);

% calculate sum of absolute correlation with empirical probality of class
p = ones(n, 1)/n; %probabilities for each classs you expect to observe in the real world

r1 =  sum(p.*abs(r));% relevance used for QPFS

F = abs(r1);
Q = abs(r2);
meanF = (1/M)*sum(F);
meanQ = (1/M^2)*sum(sum(Q));
alpha=meanQ/(meanQ + meanF);

H = (1-alpha)*Q;
f = -alpha*F;
A=[];B=[];
% sum of xi = 1
Aeq = ones(1,M);Beq = 1;

lb = zeros(1,M);
ub =[];

x = quadprog(H,f,A,B,Aeq,Beq,lb,ub);
