function tests = fun_based_unit_tests
%FUN_BASED_UNIT_TESTS   Function-based unit tests for hpl_ai_matrix.
%   This function is not meant to be called directly. To run the tests,
%   please use the script test_hpl_ai_matrix instead.
  warning('fun_based_unit_tests:invalidUsage',...
          'Please use the script ''test_hpl_ai_matrix'' tu run the tests.');
  tests = functiontests(localfunctions);
end

function test_invalid_nargin(testcase)
  verifyError(testcase,@()hpl_ai_matrix(),'hpl_ai_matrix:invalidNARGIN');
  testfun = @(x)hpl_ai_matrix(100,100,1/2,x);
  verifyError(testcase,@()testfun({}),'hpl_ai_matrix:invalidNARGIN');
  verifyError(testcase,@()testfun([]),'hpl_ai_matrix:invalidNARGIN');
  verifyError(testcase,@()testfun(1/2),'hpl_ai_matrix:invalidNARGIN');
  verifyError(testcase,@()testfun([1,2]),'hpl_ai_matrix:invalidNARGIN');
  verifyError(testcase,@()testfun('c'),'hpl_ai_matrix:invalidNARGIN');
  verifyError(testcase,@()testfun('string'),'hpl_ai_matrix:invalidNARGIN');
  verifyError(testcase,@()testfun("string"),'hpl_ai_matrix:invalidNARGIN');
  verifyError(testcase,@()testfun(@(x)(x)),'hpl_ai_matrix:invalidNARGIN');
end

function test_invalid_n(testcase)
  testfun = @(x)hpl_ai_matrix(x);
  for x = [-Inf,-1.1,-1,0,1,1i,1+1i,1.1,99,Inf,NaN]
    verifyError(testcase,@()testfun(x),'hpl_ai_matrix:invalidN');
  end
  verifyError(testcase,@()testfun({}),'hpl_ai_matrix:invalidN');
  verifyError(testcase,@()testfun([]),'hpl_ai_matrix:invalidN');
  verifyError(testcase,@()testfun([1,2]),'hpl_ai_matrix:invalidN');
  verifyError(testcase,@()testfun('c'),'hpl_ai_matrix:invalidN');
  verifyError(testcase,@()testfun('string'),'hpl_ai_matrix:invalidN');
  verifyError(testcase,@()testfun("string"),'hpl_ai_matrix:invalidN');
  verifyError(testcase,@()testfun(@(x)(x)),'hpl_ai_matrix:invalidN');
end

function test_invalid_kappa(testcase)
  testfun = @(x)hpl_ai_matrix(100,x);
  for x = [-Inf,-1.1,-1,0,0.1,0.9,1i,1+1i,Inf,NaN]
    verifyError(testcase,@()testfun(x),'hpl_ai_matrix:invalidKAPPA');
  end
  verifyError(testcase,@()testfun({}),'hpl_ai_matrix:invalidKAPPA');
  verifyError(testcase,@()testfun([]),'hpl_ai_matrix:invalidKAPPA');
  verifyError(testcase,@()testfun([1,2]),'hpl_ai_matrix:invalidKAPPA');
  verifyError(testcase,@()testfun('c'),'hpl_ai_matrix:invalidKAPPA');
  verifyError(testcase,@()testfun('string'),'hpl_ai_matrix:invalidKAPPA');
  verifyError(testcase,@()testfun("string"),'hpl_ai_matrix:invalidKAPPA');
  verifyError(testcase,@()testfun(@(x)(x)),'hpl_ai_matrix:invalidKAPPA');
end

function test_invalid_rho(testcase)
  testfun = @(x)hpl_ai_matrix(100,100,x);
  for x = [-Inf,-1.1,-1,0,1i,1+1i,1.1,Inf,NaN]
    verifyError(testcase,@()testfun(x),'hpl_ai_matrix:invalidRHO');
  end
  verifyError(testcase,@()testfun({}),'hpl_ai_matrix:invalidRHO');
  verifyError(testcase,@()testfun([]),'hpl_ai_matrix:invalidRHO');
  verifyError(testcase,@()testfun([1,2]),'hpl_ai_matrix:invalidRHO');
  verifyError(testcase,@()testfun('c'),'hpl_ai_matrix:invalidRHO');
  verifyError(testcase,@()testfun('string'),'hpl_ai_matrix:invalidRHO');
  verifyError(testcase,@()testfun("string"),'hpl_ai_matrix:invalidRHO');
  verifyError(testcase,@()testfun(@(x)(x)),'hpl_ai_matrix:invalidRHO');
end

function test_invalid_alpha(testcase)
  testfun = @(x)hpl_ai_matrix(100,[],[],x,2);
  for x = [-Inf,-1.1,-1,0,1i,1+1i,1.1,Inf,NaN]
    verifyError(testcase,@()testfun(x),'hpl_ai_matrix:invalidALPHA');
  end
  verifyError(testcase,@()testfun({}),'hpl_ai_matrix:invalidALPHA');
  verifyError(testcase,@()testfun([]),'hpl_ai_matrix:invalidALPHA');
  verifyError(testcase,@()testfun([1,2]),'hpl_ai_matrix:invalidALPHA');
  verifyError(testcase,@()testfun('c'),'hpl_ai_matrix:invalidALPHA');
  verifyError(testcase,@()testfun('string'),'hpl_ai_matrix:invalidALPHA');
  verifyError(testcase,@()testfun("string"),'hpl_ai_matrix:invalidALPHA');
  verifyError(testcase,@()testfun(@(x)(x)),'hpl_ai_matrix:invalidALPHA');
end

function test_invalid_beta(testcase)
  testfun = @(x)hpl_ai_matrix(100,[],[],1,x);
  for x = [-Inf,-1.1,-1,0,1/2,1i,1+1i,Inf,NaN]
    verifyError(testcase,@()testfun(x),'hpl_ai_matrix:invalidBETA');
  end
  verifyError(testcase,@()testfun({}),'hpl_ai_matrix:invalidBETA');
  verifyError(testcase,@()testfun([]),'hpl_ai_matrix:invalidBETA');
  verifyError(testcase,@()testfun([1,2]),'hpl_ai_matrix:invalidBETA');
  verifyError(testcase,@()testfun('c'),'hpl_ai_matrix:invalidBETA');
  verifyError(testcase,@()testfun('string'),'hpl_ai_matrix:invalidBETA');
  verifyError(testcase,@()testfun("string"),'hpl_ai_matrix:invalidBETA');
  verifyError(testcase,@()testfun(@(x)(x)),'hpl_ai_matrix:invalidBETA');
end

function test_warning_ignored_arguments(testcase)
  testfun = @(x,y)hpl_ai_matrix(100,x,y,1/2,1);
  verifyWarning(testcase,@()testfun([],1),'hpl_ai_matrix:ignoredArguments');
  verifyWarning(testcase,@()testfun(1,[]),'hpl_ai_matrix:ignoredArguments');
  verifyWarning(testcase,@()testfun(1,1),'hpl_ai_matrix:ignoredArguments');
end

function test_default_arguments(testcase)
  verifyEqual(testcase,hpl_ai_matrix(1000),hpl_ai_matrix(1000,1000,1/2));
end

function test_conditioning_pivoting(testcase)
  m = 0;
  tol = 1e-3;
  for n = [100 500 1000]
    for e = [2:14]

      kappa = 10^e;

      [A,alpha,beta] = hpl_ai_matrix(n,kappa);
      verifyLessThan(testcase,abs(cond(A,inf)-kappa)/abs(kappa),tol);
      [L,U,p] = lu(A,'vector');
      verifyEqual(testcase,p,1:n);

      % assert_eq_tol(cond(A,inf),kappa)
      % assert_eq(p,1:n)

      for rho = [1/10,1/2,4/5]

        [A,alpha,beta] = hpl_ai_matrix(n,kappa,rho);
        verifyLessThan(testcase,abs(cond(A,inf)-kappa)/abs(kappa),tol);
        [L,U,p] = lu(A,'vector');
        verifyEqual(testcase,p,1:n);

        % assert(abs(cond(A,inf)-kappa)/abs(kappa) < 1e-3);
        % assert(p,1:n)

      end
    end
  end
end

% %%%%%%%%%%%%%%%%%%%%%%%
% function assert_eq_tol(a,b)
%   m = m + 1;
%   tol = 1e-3;
%   r = abs(a-b)/abs(a);
%   if r > tol
%     error('Failure: rel difference %9.2e',r)
%   end
%   fprintf('Test %g succeeded.\n',m)
% end

% %%%%%%%%%%%%%%%%%%%%%%%
% function assert_eq(a,b)
%   m = m + 1;
%   if ~isequal(a,b)
%     A,L,U,p
%     error('Failure.')
%   end
%   fprintf('Test %g succeeded.\n',m)
% end
