"""Reaction tests

These tests specify the behavior of Reaction implementations.
"""
import unittest
import stocal


class ReactionSpecification(object) :
	def test_equality(self) :
		"""Reactions are equal if their propensity functions are identical"""
		self.assertEqual(self.Transition_A(), self.Transition_A())
		self.assertEqual(self.Transition_A(), self.Transition_B())
		self.assertNotEqual(self.Transition_A(), self.Transition_C())

	def test_propensity(self) :
		"""Inapplicable reactions have zero propensity"""
		transition = self.Transition_A()
		self.assertEqual(transition.propensity({}), 0.)

	def test_next_occurrence(self) :
		"""Inapplicable reactions occur in infinite time"""
		transition = self.Transition_A()
		self.assertEqual(transition.next_occurrence(0., {}), float('inf'))


class TestReaction(unittest.TestCase, ReactionSpecification) :
	class Transition_A(stocal.Reaction) :
		def __init__(self) :
			stocal.Reaction.__init__(self, {'a':1}, {'z':1})
		def propensity(self, state) :
			return float(state.get('a',0))

	class Transition_B(Transition_A) :
		pass

	class Transition_C(stocal.Reaction) :
		def __init__(self) :
			stocal.Reaction.__init__(self, {'a':1}, {'z':1})
		def propensity(self, state) :
			return float(state.get('a',0))


class TestMassAction(unittest.TestCase, ReactionSpecification) :
	class Transition_A(stocal.MassAction) :
		def __init__(self) :
			stocal.MassAction.__init__(self, {'a':1}, {'z':1}, 1.)

	def test_equality(self) :
		"""propensity semantics not enforced for MassAction"""
		pass


if __name__ == '__main__' :
	unittest.main()
