"""Test for stocal.examples

Most tests simply check that the example can be run without error.
"""
import unittest
import sys
from stocal.tests.test_transitions import TestReactionRule, TestMassAction

from stocal.examples.pre2017 import DegradationRule
from stocal.examples.pre2017 import LigationRule
from stocal.examples.pre2017 import AutoCatalysisRule


class TestBrusselator(unittest.TestCase):
    """Test examples.brusselator"""
    def setUp(self):
        self.stdout = sys.stdout
        sys.stdout = open('/dev/null', 'w')

    def tearDown(self):
        sys.stdout = self.stdout

    def test_example(self):
        """test running the module"""
        import stocal.examples.brusselator


class TestEvents(unittest.TestCase):
    """Test examples.events"""
    def test_example(self):
        """test process instantiation"""
        from stocal.examples.events import process
        for _ in process.trajectory({}, steps=100):
            pass


class TestPre2017(unittest.TestCase):
    """Test examples.pre2017"""
    def test_example(self):
        """test process instantiation"""
        from stocal.examples.pre2017 import process
        for _ in process.trajectory({}, steps=100):
            pass


class TestPre2017Rule(TestReactionRule):
    """Base class for rules used in pre2017"""
    def setUp(self):
        self.rule = self.Rule()

    def num_novel_reactions(self, *reactants):
        """Number of transitions returned by novel_reactions"""
        return sum(1 for _ in self.rule.novel_reactions(*reactants))


class TestPre2017DegradationRule(TestPre2017Rule):
    """Tests pre2017.DegradationRule"""
    Rule = DegradationRule

    def test_novel_reactions_number_of_reactions(self):
        """Assert correct number of inferred transitions"""
        for i in range(1, 11):
            reactant = i*'a'
            self.assertEqual(self.num_novel_reactions(reactant), (i-1))

    def test_infer_transitions_length_of_products(self):
        """Assert correct lengths of transition products"""
        for i in range(1, 11):
            reactant = i*'a'
            for trans in self.rule.novel_reactions(reactant):
                l = sum(n*len(p) for p, n in trans.products.items())
                self.assertEqual(l, i)


class TestPre2017Degradation(TestMassAction):
    """Tests pre2017.Degradation.Transition"""
    class Transition(DegradationRule.Transition):
        """Provides default constant for Reaction tests"""
        def __init__(self, reactants, products, c=1.):
            super(TestPre2017Degradation.Transition, self).__init__(
                reactants, products, c)


class TestPre2017LigationRule(TestPre2017Rule):
    """Tests for pre2017.LigationRule"""
    Rule = LigationRule

    def test_novel_reactions_number_of_reactions(self):
        """Assert correct number of inferred transitions"""
        self.assertEqual(self.num_novel_reactions('a', 'b'), 2)
        self.assertEqual(self.num_novel_reactions('a', 'a'), 2)


class TestPre2017Ligation(TestMassAction):
    """Tests for pre2017.Ligation.Transition"""
    class Transition(LigationRule.Transition):
        """Provides default constant for Reaction tests"""
        def __init__(self, reactants, products, c=1.):
            super(TestPre2017Ligation.Transition, self).__init__(
                reactants, products, c)


class TestPre2017AutoCatalysisRule(TestPre2017Rule):
    """Tests for pre2017.AutoCatalysisRule"""
    Rule = AutoCatalysisRule

    def test_novel_reactions_number_of_reactions(self):
        """Assert correct number of inferred transitions"""
        self.assertEqual(self.num_novel_reactions('a', 'b', 'aa'), 0)
        self.assertEqual(self.num_novel_reactions('a', 'b', 'ab'), 1)
        self.assertEqual(self.num_novel_reactions('a', 'b', 'ba'), 1)
        self.assertEqual(self.num_novel_reactions('a', 'a', 'aa'), 2)
        self.assertEqual(self.num_novel_reactions('ab', 'abab', 'ababab'), 2)


class TestPre2017AutoCatalysis(TestMassAction):
    """Tests for pre2017.AutoCatalysis.Transition"""
    class Transition(AutoCatalysisRule.Transition):
        """Provides default constant for Reaction tests"""
        def __init__(self, reactants, products, c=1.):
            super(TestPre2017AutoCatalysis.Transition, self).__init__(
                reactants, products, c)


if __name__ == '__main__':
    unittest.main()
