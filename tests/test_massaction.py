"""MassAction tests

These tests specify the behavior of the MassAction class.
"""
import unittest
import stocal


class MassActionSpecification(unittest.TestCase) :
	class Transition(stocal.MassAction) :
		def __init__(self, c) :
			stocal.MassAction.__init__(self, {'a':1}, {'z':1}, c)

	def test_initialization(self) :
		"""rate constant must be non-negative"""
		with self.assertRaises(ValueError) :
			self.Transition(c=-1.)

	def test_equality(self) :
		"""MassActions with different constants are considered different."""
		self.assertEqual(self.Transition(1.), self.Transition(1.))
		self.assertNotEqual(self.Transition(1.), self.Transition(1.))

	def test_equal_propensities(self) :
		"""Assert consistent propensities
		
		This test compares the propensities of a single reaction
		2*A --> 0 in a system with n molecules to a reaction system
		A+B --> 0, 2*A --> 0, 2*B -->0 in a system with n_a+n_b==n
		molecules. If the physical properties of the molecules are
		equal, the total propensity of the reactions has to be the same.
		"""
		n_a = 12; n_b=8
		c = 1.2
		a_1 = stocal.MassAction({'a':2}, {}, c).propensity({'a':n_a+n_b})
		a_2 = sum(r.propensity({'a':n_a, 'b':n_b}) for r in [
			stocal.MassAction({'a':2}, {}, c),
			stocal.MassAction({'b':2}, {}, c),
			stocal.MassAction({'a':1, 'b':1}, {}, c),
		])
		self.assertAlmostEqual(a_1, a_2)

	def test_propensities_from_equivalent_definitions(self) :
		"""Assert propensity is identical for equivalent definitions
		
		No matter whether reactants are provided as dicts or sequences,
		the propensity of the reaction must not change.
		"""
		c = 1.2
		n = 20
		a_1 = stocal.MassAction({'a': 2}, {}, c).propensity({'a': n})
		a_2 = stocal.MassAction(['a', 'a'], {}, c).propensity({'a': n})
		self.assertAlmostEqual(a_1, a_2)

if __name__ == '__main__' :
	unittest.main()
