"""Tests for the stocal package

The behaviour of stocal is specified via tests. Each stocal module
has an associated test module and each stocal class an associated
TestCase. The test suites provide AbstractTestCases for abstract base
classes. To test implementations of abstract classes, define a
TestCase that inherits from the abstract class and overload the
abstract class by the implementation class to be tested. For example,

>>> class TestMassAction(TestReaction) :
...     Transition = MassAction

would test the Reaction (and Transition) interface for the MassAction
class.
"""
