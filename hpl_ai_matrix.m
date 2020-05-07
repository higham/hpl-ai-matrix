function [A,alpha,beta] = hpl_ai_matrix(n,kappa,rho,alpha,beta)
%HPL_AI_MATRIX   Matrix for HPL-AI Benchmark.
%   A = HPL_AI_MATRIX(N,KAPPA) is an N-by-N matrix for which which
%   COND(A,inf) = KAPPA > 1 and LU factorization with partial pivoting
%   requires no row interchanges. The matrices are generated using the
%   strategy discussed in [1]. The default value of KAPPA is 1000.
%   This function is not intended for small N, and an error is returned
%   for N < 100.
%
%   [A,ALPHA,BETA] = HPL_AI_MATRIX(N,KAPPA,RHO) specifies the value of RHO
%   satisfying 0 < RHO <= 1 (default 1/2) such that the parameters ALPHA
%   and BETA defining the matrix satisfy ALPHA = RHO*BETA.
%
%   A = HPL_AI_MATRIX(N,KAPPA,RHO,ALPHA,BETA) uses the specified values
%   of ALPHA and BETA. Both parameters must be specified, and their
%   values must satisfy the two conditions 0 < ALPHA <= BETA and
%   ALPHA <= 1. KAPPA and RHO are ignored in this case, and a warning is
%   emitted if KAPPA or RHO are not empty.
%
%   Reference:
%   [1] M. Fasi and N. J. Higham. Matrices with tunable infinity-norm
%       condition number and no need for pivoting in LU factorization.
%       MIMS EPrint 2020.xx, Manchester Institute for Mathematical
%       Sciences, The University of Manchester, UK, May 2020.

  if nargin < 1
    error('hpl_ai_matrix:invalidNARGIN',...
          'At least one argument must be specified.');
  end
  if length(n) ~= 1 || ~isnumeric(n) || ~isfinite(n) ||...
        ~isreal(n) || round(n) ~= n || n < 100
      error('hpl_ai_matrix:invalidN',...
            'N must be a positive integer greater than or equal to 100.');
          end
  if nargin < 2, kappa = 1e3; end
  if (nargin < 4 && isempty(kappa)) ||...
        (nargin == 2 && ...
         (length(kappa) > 1 || ~isreal(kappa) || ~isnumeric(kappa) ||...
          ~isfinite(kappa) || kappa < 1))
    error('hpl_ai_matrix:invalidKAPPA', ...
          'KAPPA must be a finite real number greater than 1.');
  end
  if nargin < 3, rho = 1/2; end
  if (isempty(rho) && nargin < 4) ||...
        (nargin == 3 && ...
         (length(rho) > 1 || ~isreal(rho) || ~isnumeric(rho) ||...
          ~isfinite(rho) || rho <= 0 || rho > 1))
    error('hpl_ai_matrix:invalidRHO', ...
          'RHO must be a positive real less than or equal to 1.');
  end
  if nargin == 4
    error('hpl_ai_matrix:invalidNARGIN',...
          'ALPHA and BETA must both be specified.');
  end
  if nargin == 5
    if (~isempty(kappa) || ~isempty(rho))
      warning('hpl_ai_matrix:ignoredArguments',...
              'Nonempty arguments KAPPA and RHO will be ignored.');
    end
    if length(alpha) ~= 1 ||...
          ~isnumeric(alpha) || ~isfinite(alpha) || ~isreal(alpha) ||...
          alpha <= 0 || alpha > 1
      error('hpl_ai_matrix:invalidALPHA',...
            'ALHA must be a postive real smaller than 1.');
    end
    if length(beta) ~= 1 ||...
          ~isnumeric(beta) || ~isfinite(beta) || ~isreal(beta) ||...
          beta < alpha
      error('hpl_ai_matrix:invalidBETA',...
            'BETA must be a positive real no smaller than ALPHA.');
    end
  end

  if nargin <= 3

    % Compute alpha and beta to give cond(A,inf) = kappa.
    left = eps; left_val = fhpl(n,rho*left,left) - kappa;
    if left_val >= 0, error('Left endpoint should give negative f value.'), end
    right = 1/rho;
    k = 1;
    while 1
      right_val = fhpl(n,rho*right,right) - kappa;
      if isfinite(right_val) && right_val > 0, break, end
      %       fprintf('F at right endpoint, right = %9.2e, is %9.2e.\n', right, right_val)
      right = right/2;
      k = k + 1;
      if k == 100, error('Iteration cound (100) exceeded.'), end
    end
    %   fprintf('F values for [%g,%g] are [%g, %g]\n', left,right,left_val,right_val)
    beta = fzero(@(b)fhpl(n,rho*b,b)-kappa, [left right]);
    alpha = rho*beta;
    while alpha > 1
      fprintf('Initial alpha = %9.2e exceeds 1 so recomputing.\n', alpha)
      right = right/2;
      beta = fzero(@(b)fhpl(n,rho*b,b)-kappa, [eps right]);
      alpha = rho*beta;
    end

  end

  % This is effectively what we compute, but we construct it in O(N^2) flops.
  % L = gallery('triw',n,-a)';
  % U = gallery('triw',n,-b);
  % A = L*U

  % Vectorized formation.
  k = 0:n-1;
  A = diag(1 + alpha*beta*k);
  A = A + triu(-beta*ones(1,n) + alpha*beta*k',1);  % Implicit expansion.
  A = A + tril(-alpha*ones(n,1) + alpha*beta*k,-1); % Implicit expansion.

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = fhpl(n,alpha,beta)
%FHPL   Value of cond(A,inf) for matrix A(n,a,b).

  if isempty(alpha), alpha = beta/2; end
  m = length(beta);

  for p = 1:m  % Vectorize the function, to allow for use with FZERO.
    a = alpha(p); b = beta(p);

    lambda_1 = 1 + (n-1)*b;

    idash = min(floor(1/a),n);
    k = min(floor((1+b)/b),n-1);
    lambda_idash = 1 + (2*k-idash+1)*a + (n-idash)*b + ...
        (-k^2+k+3*idash*(idash-1)/2 - n*idash+n)*a*b;

    lambda_n = 1 + (2*k-n+1)*a + (-k^2+k+ n*(n-1)/2)*a*b;
    nA_est = max([lambda_1, lambda_idash, lambda_n]);

    r = (1+a)*(1+b);

    i = 1;
    delta1 = (1+a)^i*( (1 + a)^(-1) + b*(1 - r^(n-i))/(1 - r) );

    i = n;
    deltan = (1+a)^i*( (1 + a)^(-1) + b*(1 - r^(n-i))/(1 - r) );

    ninvA_est = max(delta1,deltan);

    f(p) = nA_est*ninvA_est;
    if isinf(f(p)), f(p) = realmax; end

  end
end