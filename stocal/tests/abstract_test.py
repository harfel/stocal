"""Provide tests for abstract base classes

To test implementations of abstract classes, define a
TestCase that inherits from the abstract class and overload the
abstract class by the implementation class to be tested. For example,

class TestTransition(AbstractTestCase("Transition", stocal.Transition)):
    def test_interface(self):
        "only called when self.Transition is a concrete class."
        pass

class TestMassAction(TestTransition):
    Transition = MassAction

would test the Transition interface for the MassAction class.
"""
import unittest
import inspect

def AbstractTestCase(name, cls):
    """Support tests for abstract base classes.

    To be used as base class when defining test cases for abstract
    class implementations. cls will be bound to the attribute `name`
    in the returned base class. This allows tests in the subclass to
    access the abstract base class or its concretization.
    """
    class BaseTestCase(unittest.TestCase):
        """TestCase that is skipped if the tested class is abstract."""
        def run(self, *args, **opts):
            """Run the test case only for non-abstract test classes."""
            if inspect.isabstract(getattr(self, name)):
                return
            else:
                return super(BaseTestCase, self).run(*args, **opts)
    setattr(BaseTestCase, name, cls)
    return BaseTestCase
