%% Exercise for statistical testing
%% Load example data

X_org = load('NBLexpression.dat');

z = X_org(1,:);
X = X_org(2:end,:);
size(X)

%% Exercise 1. Try t-test
% Note: ttest2 is a function in matlab/statistics-toolbox See 'help ttest2' for a detailed usage.

geneIdx=1;
X0= X(geneIdx,z==0);
X1= X(geneIdx,z==1);
[h, p, ci, stats] = ttest2( X0, X1, 'alpha', 0.01, 'dim',2 )


%% Exercise 2. See p-value histgram
% Calculate t-test statistics for 1000 genes simultaneously.

geneIdx = 1:1000;

X0= X(geneIdx,z==0);
X1= X(geneIdx,z==1);
[h, p, ci, stats] = ttest2( X0, X1, 'alpha', 0.01, 'dim',2 );

w = 0.01;
bin = w/2:w:1;
subplot(2,1,1)
hist( p, bin )
xlabel( 'P-value' )
ylabel( 'Frequency' )
subplot(2,1,2)
hist( stats.tstat,50 )
xlabel( 't-stat')
ylabel( 'Frequency' )

%% Exercise 3. Find the most significant genes
[dum, idx] = sort( p );

subplot(2,2,1)
plot( 1:1000, p(idx) )
xlabel( 'Sorted gene index' )
ylabel( 'P-value' )

gidx1 = idx(1); gidx2 = idx(2);
subplot(2,2,2)
plot( X(gidx1,z==0),X(gidx2,z==0), 'b.', ...
    X(gidx1,z==1),X(gidx2,z==1), 'rx', 'MarkerSize', 10)
xlabel('the 1st gene')
ylabel('the 2nd gene')

gidx_worst1 = idx(1000); gidx_worst2 = idx(999);
subplot(2,2,4)
plot( X(gidx_worst1,z==0),X(gidx_worst2,z==0), 'b.', ...
    X(gidx_worst1,z==1),X(gidx_worst2,z==1), 'rx', 'MarkerSize', 10)
xlabel('the worst gene')
ylabel('the 2nd worst gene')


%% Excercise 4. T-test with sub-samples
% Run the following procedure with
% (1) various random seed (i.e. repeat several times)
% (2) various geneIdx
%    such as geneIdx = gidx2, gidx_worst1, gidx_worst2

N_all = 136   % number of patients
N_sub = 20    % number of sub-sampled patients

dum = randperm(N_all);
sidx = dum( 1:N_sub );
z_sub = z(sidx);
X_sub = X(:,sidx);

geneIdx = gidx1
X0= X_sub(geneIdx,z_sub==0);
X1= X_sub(geneIdx,z_sub==1);
[h4, p4, ci4, stats4] = ttest2( X0, X1, 'alpha', 0.01, 'dim',2 )



%% Exercise 5. t-test of simulated null samples
% Exercise 2 for non-informative noise.

Xnull = randn( size(X) );

X0= Xnull(geneIdx,z==0);
X1= Xnull(geneIdx,z==1);
[h, p, ci, stats] = ttest2( X0, X1, 'alpha', 0.01, 'dim',2 );

w = 0.01;
bin = w/2:w:1;
subplot(2,1,1)
hist( p, bin )
xlabel( 'P-value' )
ylabel( 'Frequency' )
subplot(2,1,2)
hist( stats.tstat,50 )
xlabel( 't-stat')
ylabel( 'Frequency' )

%% Exercise 6. q-value analysis
% Use qvalue.m provided by Oba.

[q, pi0] = qvalue( p );
% pi0 is an estimated ratio of true null hypotheses
pi0
disp( sprintf( 'Num genes (q<0.1) : %d', sum(q<0.1) ) )
disp( sprintf( 'Num genes (q<0.2) : %d', sum(q<0.2) ) )
figure
plot( p, q, '.' )
xlabel( 'P-value' )
ylabel( 'Q-value' )

