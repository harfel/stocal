"""Process tests

These tests specify the behavior of Process and TrajectorySamplers
"""
import unittest
import stocal


class TestProcess(unittest.TestCase) :
	def test_initialization(self) :
		"""Processes can be initialized with and without reactions/rules"""
		stocal.Process([])
		stocal.Process([], [])
		stocal.Process([], rules=[])
		stocal.Process(rules=[])
		stocal.Process(transitions=[])

	def test_trajectory(self) :
		"""Process.trajectory can be called with optional arguments"""
		proc = stocal.Process([])
		self.assert_(isinstance(proc.trajectory({}), stocal.TrajectorySampler))
		self.assert_(isinstance(proc.trajectory({}, 1.), stocal.TrajectorySampler))
		self.assert_(isinstance(proc.trajectory({}, 1., 2.), stocal.TrajectorySampler))
		self.assert_(isinstance(proc.trajectory({}, t=1.), stocal.TrajectorySampler))
		self.assert_(isinstance(proc.trajectory({}, tmax=2.), stocal.TrajectorySampler))
		self.assert_(isinstance(proc.trajectory({}, steps=10), stocal.TrajectorySampler))

	def test_determinstic_trajectory(self) :
		"""Partly deterministic processes return an appropriate sapmler"""
		proc = stocal.Process([stocal.Event({}, {'a':1},1.)])
		proc.trajectory({})


class TestStaticProcess(object) :
	def test_empty_process(self) :
		"""If no transition is applicable, the sampler stops"""
		process = stocal.Process([])
		sampler = self.Sampler(process, {})
		try :
			iter(sampler).next()
			self.fail("StopIteration not raised")
		except StopIteration :
			pass
		# XXX unclear what sampler.time should be: unchanged or float('inf') ?

	def test_reaction_identity(self) :
		"""Samplers yield one of the provided transitions"""
		reaction = stocal.MassAction({},{'z':1}, 1.)
		process = stocal.Process([reaction])
		sampler = self.Sampler(process, {})
		self.assertEqual(iter(sampler).next(), reaction)
		
	def test_time_advances(self) :
		"""When yielding a transiton, time advances"""
		reaction = stocal.MassAction({},{'z':1}, 1.)
		process = stocal.Process([reaction])
		sampler = self.Sampler(process, {})
		iter(sampler).next()
		self.assertNotEqual(sampler.time, 0.)

	def test_step_advances(self) :
		"""When yielding a transiton, step advances"""
		reaction = stocal.MassAction({},{'z':1}, 1.)
		process = stocal.Process([reaction])
		sampler = self.Sampler(process, {})
		iter(sampler).next()
		self.assertEqual(sampler.step, 1)

	def test_tmax(self) :
		"""Sampler eventually advances to tmax"""
		reaction = stocal.MassAction({},{'z':1}, 1.)
		process = stocal.Process([reaction])
		sampler = self.Sampler(process, {}, tmax=5.)
		for transition in sampler : pass
		self.assertEqual(sampler.time, 5.)

	def test_steps(self) :
		"""Sampler.step advances with iteration"""
		reaction = stocal.MassAction({},{'z':1}, 1.)
		process = stocal.Process([reaction])
		sampler = self.Sampler(process, {}, steps=5)
		for transition in sampler : pass
		self.assertEqual(sampler.step, 5)

	def test_update_state(self) :
		"""Modifying the state can enable new transitions"""
		reaction = stocal.MassAction({'a':1},{}, 1.)
		process = stocal.Process([reaction])
		sampler = self.Sampler(process, {})
		for trans in sampler : pass
		sampler.update_state({'a':5})
		self.assertEqual(iter(sampler).next(), reaction)


class TestRuleProcess(object) :
	class Rule(stocal.Rule) :
		"""Degradation rule"""
		Transition = stocal.MassAction

		def infer_transitions(self, new_species, state) :
			for species in new_species :
				yield self.Transition({species:1}, {}, 1.)

	def test_dynamic_reactions(self) :
		"""One transition can generate another"""
		reaction = stocal.MassAction({'a':1}, {'b':1}, 1.)
		process = stocal.Process([reaction], [self.Rule()])
		sampler = self.Sampler(process, {'a':1})
		traj = iter(sampler)
		self.assertIs(traj.next(), reaction)
		try : traj.next()
		except StopIteration : self.fail("Transition not generated")

	def test_dynamic_state_update(self) :
		"""Modifying the state can generate new transitions through rules"""
		process = stocal.Process(rules=[self.Rule()])
		sampler = self.Sampler(process, {})
		for trans in sampler : pass
		sampler.update_state({'a':5})
		try : iter(sampler).next()
		except StopIteration : self.fail("Transition not generated")


class TestStaticDeterministicProcess(object) :
	def test_instantiation(self) :
		process = stocal.Process([
			stocal.Event({},{'a':1}, .33, 1.),
			stocal.Event({'a':1},{'b':1}, .67, 1.),
		])
		sampler = self.Sampler(process, {})


class TestDirectMethod(unittest.TestCase, TestStaticProcess, TestRuleProcess) :
	Sampler = stocal.DirectMethod


class TestFirstReactionMethod(
	unittest.TestCase,
	TestStaticProcess,
	TestRuleProcess,
	TestStaticDeterministicProcess
) :
	Sampler = stocal.FirstReactionMethod


if __name__ == '__main__' :
	unittest.main()
