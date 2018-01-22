"""Transition tests

These tests specify the behavior of Transition implementations.
"""
import unittest
import stocal


class TransitionSpecification(object) :
	def test_positive_initialization(self) :
		"""Reactants and products must be positive"""
		self.Transition({'a':1}, {'z':1})
		with self.assertRaises(ValueError) :
			self.Transition({'a':0}, {'z':1})
		with self.assertRaises(ValueError) :
			self.Transition({'a':-1}, {'z':1})
		with self.assertRaises(ValueError) :
			self.Transition({'a':1}, {'z':0})
		with self.assertRaises(ValueError) :
			self.Transition({'a':1}, {'z':-1})

	def test_nonempty_initialization(self) :
		"""Either reactants or products must be specified"""
		self.Transition({'a':1}, {})
		self.Transition({}, {'z':1})
		with self.assertRaises(ValueError) :
			self.Transition({}, {})

	def test_sequence_initialization(self) :
		"""Transitions can be intialized with dicts or sequences"""
		self.assertEqual(self.Transition(['a'],['z']), self.Transition({'a':1},{'z':1}))
		self.assertEqual(self.Transition(['a', 'a'],['z']), self.Transition({'a':2},{'z':1}))
	
	def test_equality(self) :
		"""Transitions are equal if their reactants and products are equal"""
		trans = self.Transition({'a':1},{'z':1})
		self.assertEqual(trans, self.Transition({'a':1},{'z':1}))
		self.assertNotEqual(trans, self.Transition({'b':1},{'z':1}))
		self.assertNotEqual(trans, self.Transition({'a':1},{'y':1}))

	def test_true_reactants(self) :
		"""true_reactants and true_products exclude catalysts"""
		trans1 = self.Transition({'a':1},{'a':2})
		trans2 = self.Transition({},{'a':1})
		self.assertNotEqual(trans1.true_reactants, trans1.reactants)
		self.assertNotEqual(trans1.true_products, trans1.products)
		self.assertEqual(trans1.true_reactants, trans2.true_reactants)
		self.assertEqual(trans1.true_products, trans2.true_products)

	def test_hash(self) :
		"""Equal transitions must have equal hash values"""
		trans_1 = self.Transition({'a':1},{'z':1})
		trans_2 = self.Transition({'a':1},{'z':1})
		self.assertEqual(hash(trans_1), hash(trans_2))


class TestReaction(unittest.TestCase, TransitionSpecification) :
	class Transition(stocal.Reaction) :
		def propensity(self, state) :
			return 0.

class TestMassAction(unittest.TestCase, TransitionSpecification) :
	class Transition(stocal.MassAction) :
		def __init__(self, reactants, products) :
			stocal.MassAction.__init__(self, reactants, products, 1.)

class TestEvent(unittest.TestCase, TransitionSpecification) :
	class Transition(stocal.Event) :
		def __init__(self, reactants, products) :
			stocal.Event.__init__(self, reactants, products, 1.)


if __name__ == '__main__' :
	unittest.main()
