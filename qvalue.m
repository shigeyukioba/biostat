function [q,pi0] = qvalue( p, lambda )
%  q = qvalue( p [, lambda] )
% calculate q-values from p-value distribution
% input:
%   p: column vector of p-values
%   lambda: the value of the tuning parameter to estimate pi0
%
% output:
%   q: column vector of q-values
% pi0: estimated proportion of true null hypotheses

if nargin == 1
  %  lambda = 0.8;
  lambda = [0:0.05:0.95];
end

m = length(p);

% estimate pi0
if length(lambda)==1
  pi0 = mean( p >= lambda )/(1-lambda);
else
  for i=1:length(lambda)
    pi00(i) = mean( p >= lambda(i) )/(1-lambda(i));
  end
  switch 'bootstrap'
   case 'smoother'
    disp('not implemented!')
   case 'bootstrap'
    minpi0 = min(pi00);
    mse = pi00*0;pi00boot = pi00*0;
    for bs = 1:100
      pboot = p( ceil(rand(size(p))*length(p)) );
      for i=1:length(lambda)
	pi00boot(i) = mean( pboot >= lambda(i) )/(1-lambda(i));
      end
      mse = mse + (pi00boot-minpi0).^2;
    end
    pi0 = min( pi00(mse==min(mse)) );
  end
end

% calculate q-values

v = ranking(p);
q = min( pi0 * m * p ./ v, 1.0 );


function r = ranking(x)
n = length(x);
[dum, idx] = sort(x);
r = zeros(n,1);
r(idx) = 1:n;
