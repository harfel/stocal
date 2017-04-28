"""Reaction tests

These tests specify the behavior of Reaction and its
standard subclass MassAction.
"""
import unittest
import stocal


class TestReaction(unittest.TestCase) :
	class R(stocal.Reaction) :
		def propensity(state) :
			return 1.

	class S(R) :
		pass

	class T(stocal.Reaction) :
		def propensity(state) :
			return 2.
		
	def test_initialization(self) :
		"""Init specification
		
		Reactants and products must be non-negative
		"""
		with self.assertRaises(ValueError) :
			self.R({'a':0}, {'z':1})
		with self.assertRaises(ValueError) :
			self.R({'a':-1}, {'z':1})
		with self.assertRaises(ValueError) :
			self.R({'a':1}, {'z':0})
		with self.assertRaises(ValueError) :
			self.R({'a':1}, {'z':-1})

	def test_zero_order_initialization(self) :
		"""Empty reactans or products
		
		Reactions can be initialized with empty reactant or
		product sets (but not both of them)."""
		try :
			self.R({}, {'z':1})
		except :
			self.fail("empty reactant set not supported.")
		try :
			self.R({'a':1}, {})
		except :
			self.fail("empty product set not supported.")
		with self.assertRaises(ValueError) :
			self.R({}, {})

	def test_equality(self) :
		"""Equality specification
		
		Reaction's are equal if and only if their reactants,
		products, and propensity functions are equal. 
		"""
		self.assertEqual( self.R({'a':1},{'z':1}), self.R({'a':1},{'z':1}) )
		self.assertNotEqual( self.R({'a':1},{'z':1}), self.R({'b':1},{'z':1}) )
		self.assertNotEqual( self.R({'a':1},{'z':1}), self.R({'a':1},{'y':1}) )

		self.assertEqual( self.R({'a':1},{'z':1}), self.S({'a':1},{'z':1}) )
		self.assertNotEqual( self.R({'a':1},{'z':1}), self.S({'b':1},{'z':1}) )
		self.assertNotEqual( self.R({'a':1},{'z':1}), self.S({'a':1},{'y':1}) )

		self.assertNotEqual( self.R({'a':1},{'z':1}), self.T({'a':1},{'z':1}) )
		self.assertNotEqual( self.R({'a':1},{'z':1}), self.T({'b':1},{'z':1}) )
		self.assertNotEqual( self.R({'a':1},{'z':1}), self.T({'a':1},{'y':1}) )

	def test_hashing(self) :
		"""Hash specification
		
		Reactions that are equal product the same hash value
		"""
		r1 = self.R({'a':1},{'z':1})
		r2 = self.R({'a':1},{'z':1})
		r3 = self.S({'a':1},{'z':1})
		self.assertEqual( hash(r1), hash(r2) )
		self.assertEqual( hash(r1), hash(r3) )

	@unittest.skip("not enforced, but documented.")
	def test_immutablity(self) :
		"""Immutability specification
		
		Changing reactants or products of a reaction raises an AttributeError
		"""
		r = self.R({'a':1}, {'z':1})
		with self.assertRaises(AttributeError) :
			r.reactants = {'b':1}
		with self.assertRaises(AttributeError) :
			r.products = {'y':1}
		with self.assertRaises(AttributeError) :
			r.reactants['a'] = 2
		with self.assertRaises(AttributeError) :
			r.products['z'] = 2


class TestMassAction(TestReaction) :
	class R(stocal.MassAction) :
		def __init__(self, reactants, products) :
			super(TestMassAction.R, self).__init__(reactants, products, 1.)

	class S(R) :
		pass

	class T(stocal.MassAction) :
		def __init__(self, reactants, products) :
			super(TestMassAction.T, self).__init__(reactants, products, 2.)

	def test_negative_initialization(self) :
		"""Init specification
		
		Reactants and products must be non-negative
		"""
		with self.assertRaises(ValueError) :
			stocal.MassAction({'a':1}, {'z':1}, -1.)

	def test_equality_with_reaction(self) :
		"""MassActions are not equal to Reactions"""
		self.assertNotEqual(self.R({'a':1},{'z':1}), TestReaction.R({'a':1},{'z':1}))

	def test_example_propensities(self) :
		"""Test some known propensities"""
		self.assertEqual(self.R({'a':1},{}).propensity({}), 0.)
		self.assertEqual(self.R({'a':1},{}).propensity({'b':1}), 0.)
		self.assertEqual(self.R({'a':2},{}).propensity({'a':1}), 0.)
		self.assertAlmostEqual(self.R({}, {'a':1}).propensity({}), 1.)
		
	def test_indistinguishable_propensity(self) :
		"""Test some general propensity relationships"""
		r1 = self.R({'a':2}, {'z':1}).propensity
		r2 = self.R({'a':1, 'b':1}, {'z':1}).propensity
		self.assertAlmostEqual(r1({'a':10}), 0.5*10*r2({'a':1, 'b':9}))


if __name__ == '__main__' :
	unittest.main()
