"""Process tests

These tests specify the behavior of the Process class.
"""
import unittest
import stocal


class TestProcess(unittest.TestCase) :
	def test_init(self) :
		"""Processes can be initialized with and without reactions/rules"""
		stocal.Process([])
		stocal.Process([], [])
		stocal.Process([], rules=[])
		stocal.Process(rules=[])
		stocal.Process(transitions=[])

	def test_trajectory_arguments(self) :
		"""Process.trajectory can be called with optional arguments"""
		proc = stocal.Process([])
		proc.trajectory({})
		proc.trajectory({}, 1.)
		proc.trajectory({}, 1., 2.)
		proc.trajectory({}, t=1.)
		proc.trajectory({}, tmax=2.)
		proc.trajectory({}, steps=10)

	def test_trajectory_with_events(self) :
		"""Partly deterministic processes return an appropriate sapmler"""
		proc = stocal.Process([stocal.Event({}, {'a':1},1.)])
		proc.trajectory({})


if __name__ == '__main__' :
	unittest.main()
