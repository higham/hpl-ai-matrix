% This script runs the unit tests for the hpl_ai_matrix function.
warning('off', 'fun_based_unit_tests:invalidUsage')
test_results = run(fun_based_unit_tests);
warning('on', 'fun_based_unit_tests:invalidUsage')
table(test_results)
