import unittest
import sys
import stocal
from test_transitions import TestReactionRule, TestMassAction


class TestBrusselator(unittest.TestCase) :
	def setUp(self) :
		self.stdout = sys.stdout
		sys.stdout = open('/dev/null', 'w')

	def tearDown(self) :
		sys.stdout = self.stdout	

	def test_example(self) :
		import stocal.examples.brusselator


class TestEvents(unittest.TestCase) :
	def test_example(self) :
		from stocal.examples.events import process
		for _ in process.trajectory({}, steps=100) :
			pass


from stocal.examples.pre2017 import DegradationRule
from stocal.examples.pre2017 import LigationRule
from stocal.examples.pre2017 import AutoCatalysisRule

class TestPre2017(unittest.TestCase) :
	def test_example(self) :
		from stocal.examples.pre2017 import process
		for _ in process.trajectory({}, steps=100) :
			pass


class TestPre2017Rule(TestReactionRule) :
	def setUp(self) :
		self.rule = self.Rule()

	def num_novel_reactions(self, *reactants) :
		return sum(1 for _ in self.rule.novel_reactions(*reactants))

class TestPre2017DegradationRule(TestPre2017Rule) :
	Rule = DegradationRule

	def test_novel_reactions_number_of_reactions(self) :
		for i in range(1, 11):
			reactant = i*'a'
			self.assertEqual(self.num_novel_reactions(reactant), (i-1))

	def test_infer_transitions_length_of_products(self) :
		for i in range(1, 11):
			reactant = i*'a'
			for trans in self.rule.novel_reactions(reactant) :
				l = sum(n*len(p) for p,n in trans.products.iteritems())
				self.assertEqual(l, i)


class TestPre2017Degradation(TestMassAction) :
	class Transition(DegradationRule.Transition) :
		def __init__(self, reactants, products, c=1.) :
			super(TestPre2017Degradation.Transition, self).__init__(
			      reactants, products, c)


class TestPre2017LigationRule(TestPre2017Rule) :
	Rule = LigationRule

	def test_novel_reactions_number_of_reactions(self) :
		self.assertEqual(self.num_novel_reactions('a', 'b'), 2)
		self.assertEqual(self.num_novel_reactions('a', 'a'), 2)


class TestPre2017Ligation(TestMassAction) :
	class Transition(LigationRule.Transition) :
		def __init__(self, reactants, products, c=1.) :
			super(TestPre2017Ligation.Transition, self).__init__(
			      reactants, products, c)


class TestPre2017AutoCatalysisRule(TestPre2017Rule) :
	Rule = AutoCatalysisRule

	def test_novel_reactions_number_of_reactions(self) :
		self.assertEqual(self.num_novel_reactions('a', 'b', 'aa'), 0)
		self.assertEqual(self.num_novel_reactions('a', 'b', 'ab'), 1)
		self.assertEqual(self.num_novel_reactions('a', 'b', 'ba'), 1)
		self.assertEqual(self.num_novel_reactions('a', 'a', 'aa'), 2)
		self.assertEqual(self.num_novel_reactions('ab', 'abab', 'ababab'), 2)


class TestPre2017AutoCatalysis(TestMassAction) :
	class Transition(AutoCatalysisRule.Transition) :
		def __init__(self, reactants, products, c=1.) :
			super(TestPre2017AutoCatalysis.Transition, self).__init__(
			      reactants, products, c)


if __name__ == '__main__' :
	unittest.main()
