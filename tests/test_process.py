"""Process tests

These tests specify the behavior of Process and TrajectorySamplers
"""
import unittest
import stocal


class TestStaticProcess(unittest.TestCase) :
	def test_trajectory(self) :
		"""Process.trajectory can be called with optional arguments"""
		proc = stocal.Process([])
		self.assert_(isinstance(proc.trajectory({}), stocal.TrajectorySampler))
		self.assert_(isinstance(proc.trajectory({}, 1.), stocal.TrajectorySampler))
		self.assert_(isinstance(proc.trajectory({}, 1., 2.), stocal.TrajectorySampler))
		self.assert_(isinstance(proc.trajectory({}, t=1.), stocal.TrajectorySampler))
		self.assert_(isinstance(proc.trajectory({}, tmax=2.), stocal.TrajectorySampler))
		self.assert_(isinstance(proc.trajectory({}, steps=10), stocal.TrajectorySampler))


class TestTrajectorySampler(object) :
	def test_empty_process(self) :
		"""If no transition is applicable, the sampler stops"""
		process = stocal.Process([])
		traj = process.trajectory({})
		try :
			iter(traj).next()
			self.fail("StopIteration not raised")
		except StopIteration :
			pass

	def test_reaction_identity(self) :
		"""Samplers yield one of the provided transitions"""
		reaction = stocal.MassAction({},{'z':1}, 1.)
		process = stocal.Process([reaction])
		traj = process.trajectory({})
		self.assertEqual(iter(traj).next(), reaction)
		
	def test_time_advances(self) :
		"""When yielding a transiton, time advances"""
		reaction = stocal.MassAction({},{'z':1}, 1.)
		process = stocal.Process([reaction])
		traj = process.trajectory({})
		iter(traj).next()
		self.assertNotEqual(traj.time, 0.)

	def test_step_advances(self) :
		"""When yielding a transiton, step advances"""
		reaction = stocal.MassAction({},{'z':1}, 1.)
		process = stocal.Process([reaction])
		traj = process.trajectory({})
		iter(traj).next()
		self.assertEqual(traj.step, 1)

	def test_tmax(self) :
		"""Sampler eventually advances to tmax"""
		reaction = stocal.MassAction({},{'z':1}, 1.)
		process = stocal.Process([reaction])
		traj = process.trajectory({}, tmax=5.)
		for transition in traj : pass
		self.assertEqual(traj.time, 5.)

	def test_steps(self) :
		reaction = stocal.MassAction({},{'z':1}, 1.)
		process = stocal.Process([reaction])
		traj = process.trajectory({}, steps=5)
		for transition in traj : pass
		self.assertEqual(traj.step, 5)

	def test_update_state(self) :
		"""Modifying the state can enable new transitions"""
		reaction = stocal.MassAction({'a':1},{}, 1.)
		process = stocal.Process([reaction])
		traj = process.trajectory({})
		try :
			iter(traj).next()
			self.fail("StopIteration not raised")
		except StopIteration :
			pass
		traj.update_state({'a':5})
		self.assertEqual(iter(traj).next(), reaction)


class TestDirectMethod(unittest.TestCase, TestTrajectorySampler) :
	Sampler = stocal.DirectMethod
