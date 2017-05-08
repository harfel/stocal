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

	def test_propensities(self) :
		"""Assert consistent propensities
		
		This test compares the propensities of a single reaction
		2*A --> 0 in a system with n molecules to a reaction system
		A+B --> 0, 2*A --> 0, 2*B -->0 in a system with n_a+n_b==n
		molecules.
		"""
		n_a = 12; n_b=8
		p = 1.2
		p_1 = stocal.MassAction({'a':2}, {}, p).propensity({'a':n_a+n_b})
		p_2 = sum(r.propensity({'a':n_a, 'b':n_b}) for r in [
			stocal.MassAction({'a':2}, {}, p),
			stocal.MassAction({'b':2}, {}, p),
			stocal.MassAction({'a':1, 'b':1}, {}, p),
		])
		self.assertAlmostEqual(p_1, p_2)


if __name__ == '__main__' :
	unittest.main()
