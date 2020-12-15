"""
Test to determine if the behavior exhibited by the new check_cpp_support() function matches
what is expected by the previous try/catch implementation.
"""

import unittest;

class TestCheckCPPSupport(unittest.TestCase):
    def test_check_cpp_support(self):
        from gillespy2.solvers.utilities.cpp_support_test import check_cpp_support
        self.assertEqual(check_cpp_support(), self.old_check_cpp_support())

    def old_check_cpp_support(self):
        from gillespy2.solvers.cpp.example_models import Example
        from gillespy2 import SSACSolver
        try:
            model = Example()
            results = model.run(solver=SSACSolver, cpp_support=True)
            return True
        except Exception as e:
            from gillespy2.core import log
            log.warn('Unable to use C++ optimized SSA: {0}.  The performance of ' \
            'this package can be significantly increased if you install/configure GCC on ' \
            'this machine.'.format(e))
            return False

if __name__ == '__main__':
    unittest.main()
