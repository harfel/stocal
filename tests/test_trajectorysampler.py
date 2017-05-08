"""TrajectorySampler tests

These tests specify the behavior of TrajectorySampler implementations.
"""
import unittest
import stocal


class TrajectorySamplerSpecification(object) :
	def test_init_optional_args(self) :
		"""init can be called with optional arguments"""
		proc = stocal.Process([])
		self.assertEqual(self.Sampler(proc, {}).time, 0.)
		self.assertEqual(self.Sampler(proc, {}, t=1.).time, 1.)
		self.assertEqual(self.Sampler(proc, {}, tmax=10.).tmax, 10.)
		self.assertEqual(self.Sampler(proc, {}, steps=10).steps, 10)

	def test_init_bounds(self) :
		"""accepted parameters are within range
		
		state must have only positive copy numbers
		t must be non-negative
		tmax must be non-negative
		steps must non-negative
		"""
		proc = stocal.Process([])
		with self.assertRaises(ValueError) :
			self.Sampler(proc, {}, t=-1.)
		with self.assertRaises(ValueError) :
			self.Sampler(proc, {}, tmax=-1.)
		with self.assertRaises(ValueError) :
			self.Sampler(proc, {}, steps=-1)
		with self.assertRaises(ValueError) :
			self.Sampler(proc, {'a':0})
		with self.assertRaises(ValueError) :
			self.Sampler(proc, {'a':-1})

	def test_add_transition_enables_transition(self) :
		"""transitions added with add_transition must be executable"""
		sampler = self.Sampler(stocal.Process([]), {'a':1})
		transition = stocal.MassAction({'a':1}, {}, 1.)
		sampler.add_transition(transition)
		self.assertIs(iter(sampler).next(), transition)

	def test_update_state_enables_static(self) :
		"""update_state can enable static transitions"""
		transition = stocal.MassAction({'a':1}, {}, 1.)
		sampler = self.Sampler(stocal.Process([transition]), {})
		sampler.update_state({'a':1})
		self.assertIs(iter(sampler).next(), transition)

	def test_update_state_disables_static(self) :
		"""update_state can disable static transitions"""
		transition = stocal.MassAction({'a':1}, {}, 1.)
		sampler = self.Sampler(stocal.Process([transition]), {'a':1})
		sampler.update_state({'a':0})
		with self.assertRaises(StopIteration) :
			iter(sampler).next()

	def test_update_state_enables_infered(self) :
		"""update_state can enable infered transitions"""
		class Rule(stocal.Rule) :
			Transition = stocal.MassAction
			def infer_transitions(self, new_species, state) :
				for species in new_species :
					yield self.Transition({species:1}, {}, 1.)

		sampler = self.Sampler(stocal.Process([],[Rule()]), {})
		sampler.update_state({'a':1})
		self.assert_(isinstance(iter(sampler).next(), Rule.Transition))

	def test_update_state_disables_infered(self) :
		"""update_state can disable infered transitions"""
		class Rule(stocal.Rule) :
			Transition = stocal.MassAction
			def infer_transitions(self, new_species, state) :
				for species in new_species :
					yield self.Transition({species:1}, {}, 1.)

		sampler = self.Sampler(stocal.Process([],[Rule()]), {'a':1})
		sampler.update_state({'a':0})
		with self.assertRaises(StopIteration) :
			iter(sampler).next()

	def test_iter_empty(self) :
		"""The empty process stops iteration immediately"""
		sampler = self.Sampler(stocal.Process([]), {'a':100})
		with self.assertRaises(StopIteration) :
			iter(sampler).next()
		self.assertEqual(sampler.step, 0)
		self.assertEqual(sampler.time, 0.)

	def test_iter_empty_tmax(self) :
		"""If tmax is given, sampler.time advances to it"""
		sampler = self.Sampler(stocal.Process([]), {'a':100}, tmax=100.)
		with self.assertRaises(StopIteration) :
			iter(sampler).next()
		self.assertEqual(sampler.step, 0)
		self.assertEqual(sampler.time, 100.)

	def test_iter_steps(self) :
		"""If steps is given, sampler returns this many transitions"""
		n = 50
		transition = stocal.MassAction({'a':1}, {}, 1.)
		sampler = self.Sampler(stocal.Process([transition]), {'a':100}, steps=n)
		self.assertEqual(sum(1 for _ in sampler), n)
		self.assertEqual(sampler.steps, n)

	def test_iter_zero_steps(self) :
		"""If steps equals zero, sampler stops immediately"""
		n = 0
		transition = stocal.MassAction({'a':1}, {}, 1.)
		sampler = self.Sampler(stocal.Process([transition]), {'a':100}, steps=n)
		self.assertEqual(sum(1 for _ in sampler), n)
		self.assertEqual(sampler.steps, n)

	def test_iter_steps_and_tmax(self) :
		"""if steps exceed, time should not be advanced to tmax"""
		transition = stocal.MassAction({'a':1}, {}, 1.)
		proc = stocal.Process([transition])
		sampler = self.Sampler(proc, {'a':100}, steps=0, tmax=10.)
		with self.assertRaises(StopIteration) :
			iter(sampler).next()
		self.assertEqual(sampler.step, 0)
		self.assertEqual(sampler.time, 0.)


class TestDirectMethod(unittest.TestCase, TrajectorySamplerSpecification) :
	Sampler = stocal.DirectMethod


class FirstReactionMethodAction(unittest.TestCase, TrajectorySamplerSpecification) :
	Sampler = stocal.FirstReactionMethod

	def test_iter_simultaneous_events(self) :
		"""two Events can occur at the exact same time"""
		proc = stocal.Process([
			stocal.Event({},{'a':1}, 10),
			stocal.Event({},{'b':1}, 10),
		])
		sampler = self.Sampler(proc, {})
		for trans in sampler : pass
		self.assertEqual(sampler.state, {'a':1, 'b':1})
		self.assertEqual(sampler.step, 2)
		self.assertEqual(sampler.time, 10)


if __name__ == '__main__' :
	unittest.main()
