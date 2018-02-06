"""Tests for stocal.algorithms"""
import unittest
from stocal.tests.abstract_test import AbstractTestCase
import stocal


class TestTrajectorySampler(AbstractTestCase('Sampler', stocal.algorithms.TrajectorySampler)):
    """Test Trajectory sampler interface

	This abstract test case provides tests that specify the
	TrajectorySampler interface. Every concrete TrajectorySampler should
	provide a test case derived from TestTrajectorySampler which stes
	the class attribute Sampler to the implementation class.
    """
    def test_init_optional_args(self):
        """init can be called with optional arguments"""
        proc = stocal.Process([])
        self.assertEqual(self.Sampler(proc, {}).time, 0.)
        self.assertEqual(self.Sampler(proc, {}, t=1.).time, 1.)
        self.assertEqual(self.Sampler(proc, {}, tmax=10.).tmax, 10.)
        self.assertEqual(self.Sampler(proc, {}, steps=10).steps, 10)

    def test_init_bounds(self):
        """accepted parameters are within range

        state must have only positive copy numbers
        t must be non-negative
        tmax must be non-negative
        steps must non-negative
        """
        proc = stocal.Process([])
        with self.assertRaises(ValueError):
            self.Sampler(proc, {}, t=-1.)
        with self.assertRaises(ValueError):
            self.Sampler(proc, {}, tmax=-1.)
        with self.assertRaises(ValueError):
            self.Sampler(proc, {}, steps=-1)
        with self.assertRaises(ValueError):
            self.Sampler(proc, {'a':0})
        with self.assertRaises(ValueError):
            self.Sampler(proc, {'a':-1})

    def test_add_transition_enables_transition(self):
        """transitions added with add_transition must be executable"""
        sampler = self.Sampler(stocal.Process([]), {'a':1})
        transition = stocal.MassAction({'a':1}, {}, 1.)
        sampler.add_transition(transition)
        self.assertIs(next(iter(sampler)), transition)

    def test_update_state_enables_static(self):
        """update_state can enable static transitions"""
        transition = stocal.MassAction({'a':1}, {}, 1.)
        sampler = self.Sampler(stocal.Process([transition]), {})
        sampler.update_state({'a':1})
        self.assertIs(next(iter(sampler)), transition)

    def test_update_state_disables_static(self):
        """update_state can disable static transitions"""
        transition = stocal.MassAction({'a':1}, {}, 1.)
        sampler = self.Sampler(stocal.Process([transition]), {'a':1})
        sampler.update_state({'a':0})
        with self.assertRaises(StopIteration):
            next(iter(sampler))

    def test_update_state_enables_infered(self):
        """update_state can enable infered transitions"""
        class Rule(stocal.Rule):
            """Rule that infers one reaction for any species"""
            Transition = stocal.MassAction
            def infer_transitions(self, new_species, state):
                for species in new_species:
                    yield self.Transition({species:1}, {}, 1.)

        sampler = self.Sampler(stocal.Process([], [Rule()]), {})
        sampler.update_state({'a':1})
        self.assertTrue(isinstance(next(iter(sampler)), Rule.Transition))

    def test_update_state_disables_infered(self):
        """update_state can disable infered transitions"""
        class Rule(stocal.Rule):
            """Rule that infers one reaction for any species"""
            Transition = stocal.MassAction
            def infer_transitions(self, new_species, state):
                for species in new_species:
                    yield self.Transition({species:1}, {}, 1.)

        sampler = self.Sampler(stocal.Process([], [Rule()]), {'a':1})
        sampler.update_state({'a':0})
        with self.assertRaises(StopIteration):
            next(iter(sampler))

    def test_iter_empty(self):
        """The empty process stops iteration immediately"""
        sampler = self.Sampler(stocal.Process([]), {'a':100})
        with self.assertRaises(StopIteration):
            next(iter(sampler))
        self.assertEqual(sampler.step, 0)
        self.assertEqual(sampler.time, 0.)

    def test_iter_empty_tmax(self):
        """If tmax is given, sampler.time advances to it"""
        sampler = self.Sampler(stocal.Process([]), {'a':100}, tmax=100.)
        with self.assertRaises(StopIteration):
            next(iter(sampler))
        self.assertEqual(sampler.step, 0)
        self.assertEqual(sampler.time, 100.)

    def test_iter_steps(self, steps=50):
        """If steps is given, sampler returns this many transitions"""
        transition = stocal.MassAction({'a':1}, {}, 1.)
        sampler = self.Sampler(stocal.Process([transition]), {'a':100}, steps=steps)
        self.assertEqual(sum(1 for _ in sampler), steps)
        self.assertEqual(sampler.step, steps)

    def test_iter_zero_steps(self):
        """If steps equals zero, sampler stops immediately"""
        return self.test_iter_steps(0)

    def test_iter_steps_and_tmax(self):
        """if steps exceed, time should not be advanced to tmax"""
        transition = stocal.MassAction({'a':1}, {}, 1.)
        proc = stocal.Process([transition])
        sampler = self.Sampler(proc, {'a':100}, steps=0, tmax=10.)
        with self.assertRaises(StopIteration):
            next(iter(sampler))
        self.assertEqual(sampler.step, 0)
        self.assertEqual(sampler.time, 0.)



class TestDirectMethod(TestTrajectorySampler):
    """Test stocal.algorithms.DirectMethod

    This tests the regular TrajectorySampler interface."""
    Sampler = stocal.algorithms.DirectMethod


class TestFirstReactionMethod(TestTrajectorySampler):
    """Test stocal.algorithms.FirstReactionMethod

    In addition to regular TrajectorySampler's this implementation
    has to properly work with deterministic stocal.Event's"""
    Sampler = stocal.algorithms.FirstReactionMethod

    def test_iter_simultaneous_events(self):
        """two Events can occur at the exact same time"""
        proc = stocal.Process([
            stocal.Event({}, {'a':1}, 10),
            stocal.Event({}, {'b':1}, 10),
        ])
        sampler = self.Sampler(proc, {})
        for _ in sampler:
            pass
        self.assertEqual(sampler.state, {'a':1, 'b':1})
        self.assertEqual(sampler.step, 2)
        self.assertEqual(sampler.time, 10)

    def test_exact_number_of_events(self):
        """sampler performs specified number of events"""
        proc = stocal.Process([
            stocal.Event({}, {'a':1}, 0, 1)
        ])
        sampler = self.Sampler(proc, {}, tmax=10)
        for _ in sampler:
            pass
        self.assertEqual(sampler.state, {'a':10})
        self.assertEqual(sampler.step, 10)
        self.assertAlmostEqual(sampler.time, 10)


if __name__ == '__main__':
    unittest.main()
