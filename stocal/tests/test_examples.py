"""Test for stocal.examples

Most tests simply check that the example can be run without error.
"""
import unittest
import sys
import os

from stocal.tests.test_transitions import TestReactionRule as TestTransitionRule, TestMassAction

from stocal.examples.pre2017 import DegradationRule
from stocal.examples.pre2017 import LigationRule
from stocal.examples.pre2017 import AutoCatalysisRule


class TestBrusselator(unittest.TestCase):
    """Test examples.brusselator"""
    def test_example(self):
        """test process instantiation"""
        from stocal.examples.brusselator import process
        for _ in process.sample({}, steps=100):
            pass


class TestEvents(unittest.TestCase):
    """Test examples.events"""
    def test_example(self):
        """test process instantiation"""
        from stocal.examples.events import process
        for _ in process.sample({}, steps=100):
            pass


class TestPre2017(unittest.TestCase):
    """Test examples.pre2017"""
    def test_example(self):
        """test process instantiation"""
        from stocal.examples.pre2017 import process
        for _ in process.sample({}, steps=100):
            pass


class TestPre2017Rule(TestTransitionRule):
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


class TestTypedRules(TestTransitionRule):
    from stocal.examples.typed_rules import AA_BB as Rule

    def test_infer_transitions_signature(self):
        from stocal.examples.typed_rules import AA, BB
        
        rule = self.Rule()
        transitions = list(rule.infer_transitions({AA('a'): 2},
                                                  {BB('z'): 1}))
        self.assertEqual(len(transitions), 2)

        transitions = list(rule.infer_transitions({AA('a'): 1},
                                                  {AA('a'): 1, BB('z'): 1}))
        self.assertEqual(len(transitions), 0)


class TestTemperatureCycle(unittest.TestCase):
    """Test temperature_cycle example"""
    def test_example(self):
        """test process instantiation"""
        from stocal.examples.temperature_cycle import process
        for _ in process.sample({}, steps=100):
            pass


class TestValidation(unittest.TestCase):
    """Test validation example"""
    from stocal.examples.validation import DataStore
    
    def setUp(self):
        from tempfile import mkdtemp
        self.tmpdir = mkdtemp(prefix='stocal-tmp')
        self.store = self.DataStore(self.tmpdir)

    def tearDown(self):
        from shutil import rmtree
        rmtree(self.tmpdir)

    def test_run(self):
        """Assert that validation run is executable"""
        from argparse import Namespace
        from stocal.algorithms import DirectMethod as TestMethod
        from stocal.examples.dsmts.models import DSMTS_001_01 as TestModel
        from stocal.examples.validation import run_validation

        args = Namespace(models=[TestModel], algo=[TestMethod], N=3,
                         cpu=1, store=self.store)
        run_validation(args)

    def test_report(self):
        """Assert that validation report is executable"""
        # populate the store first...
        from argparse import Namespace
        from stocal.examples.validation import report_validation

        report_name = os.path.join(self.tmpdir, 'validation.tex')
        args = Namespace(reportfile=report_name, cpu=1, store=self.store)
        self.test_run()
        report_validation(args)


if __name__ == '__main__':
    unittest.main()
