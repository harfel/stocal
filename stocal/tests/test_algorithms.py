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

    def test_add_transition_enables_transition(self):
        """transitions added with add_transition must be executable"""
        sampler = self.Sampler(stocal.Process([]), {'a':1})
        transition = stocal.MassAction({'a':1}, {}, 1.)
        sampler.add_transition(transition)
        self.assertIs(next(iter(sampler)), transition)

    def test_propose_potential_transition_empty(self):
        """Proposed (time,transition) for empty process is (inf,None)"""
        sampler = self.Sampler(stocal.Process([]), {'a':100}, tmax=100.)
        time, transition, args = sampler.propose_potential_transition()
        self.assertEqual(time, float('inf'))
        self.assertEqual(transition, None)

    def test_propose_potential_transition_seed(self):
        """Samplers initialized with the same random seed propose equal transitions"""
        process = stocal.Process([stocal.MassAction([], ['a'], 1.)])
        sampler_a = self.Sampler(process, {'a':100}, seed=10)
        sampler_b = self.Sampler(process, sampler_a.state, seed=10)
        self.assertEqual(sampler_a.propose_potential_transition(), sampler_b.propose_potential_transition())

    def test_propose_potential_transition_in_finite_time(self):
        """Proposed (time,transition) for empty process is (inf,None)"""
        process = stocal.Process([stocal.MassAction([], ['a'], 1.)])
        sampler = self.Sampler(process, {}, tmax=100.)
        time, transition, args = sampler.propose_potential_transition()
        self.assertTrue(time < float('inf'))
        self.assertNotEqual(transition, None)

    def test_perform_transition_advances_steps(self):
        transition = stocal.MassAction([], ['a'], 1.)
        process = stocal.Process([transition])
        sampler = self.Sampler(process, {}, tmax=100.)
        time, trans, args = sampler.propose_potential_transition()
        sampler.perform_transition(time, trans, *args)
        self.assertEqual(sampler.step, 1)

    def test_perform_transition_advances_time(self):
        transition = stocal.MassAction([], ['a'], 1.)
        process = stocal.Process([transition])
        sampler = self.Sampler(process, {}, tmax=100.)
        time, trans, args = sampler.propose_potential_transition()
        sampler.perform_transition(time, trans, *args)
        self.assertGreater(sampler.time, 0.)

    def test_perform_transition_changes_state(self):
        transition = stocal.MassAction([], ['a'], 1.)
        process = stocal.Process([transition])
        sampler = self.Sampler(process, {}, tmax=100.)
        time, trans, args = sampler.propose_potential_transition()
        sampler.perform_transition(time, trans, *args)
        self.assertEqual(sampler.state, {'a':1})

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

    def test_transitions_counts_mutliplicities(self):
        """Sampler.transitions should give access to all transitions."""
        proc = stocal.Process()
        sampler = self.Sampler(proc, {})
        sampler.add_transition(stocal.MassAction({}, {'a':1}, 1.))
        sampler.add_transition(stocal.MassAction({}, {'a':1}, 1.))
        sampler.add_transition(stocal.MassAction({}, {'b':1}, 1.))
        self.assertEqual(len(sampler.transitions), 3)


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

    def test_iter_includes_all_events_at_tmax(self):
        proc = stocal.Process([
            stocal.Event({}, {'a':1}, 10, 10),
            stocal.Event({}, {'b':1},  0, 10),
        ])
        sampler = self.Sampler(proc, {}, tmax=10)
        for _ in sampler:
            pass
        self.assertEqual(sampler.step, 3)
        self.assertEqual(sampler.time, 10)
        self.assertEqual(sampler.state, {'a':1, 'b':2})

    def test_exact_number_of_events(self):
        """sampler performs specified number of events"""
        proc = stocal.Process([
            stocal.Event({}, {'a':1}, 0, 1)
        ])
        sampler = self.Sampler(proc, {}, tmax=10)
        for _ in sampler:
            pass
        self.assertEqual(sampler.step, 11)
        self.assertAlmostEqual(sampler.time, 10)

    def test_fire_inferred_event(self):
        """sampler fires inferred events"""
        class Rule(stocal.ReactionRule):
            Transition = stocal.Event
    
            def novel_reactions(self, x):
                if x=='a':
                    yield self.Transition(['a'], [], 10)

        process = stocal.Process(transitions=[stocal.Event([], ['a'], 1)],
                                 rules=[Rule()])

        traj = stocal.algorithms.AndersonNRM(process, {})
        it = iter(traj)

        try:
            trans = next(it)
        except StopIteration:
            self.fail("Static event not fired.")
        self.assertEqual(traj.time, 1)
        self.assertEqual(traj.state, {'a': 1})

        try:
            trans = next(it)
        except StopIteration:
            self.fail("Infered event not fired.")
        self.assertEqual(traj.time, 10)
        self.assertEqual(traj.state, {})

    def test_do_not_apply_inapplicable_events(self):
        """assert that Event does not fire if reactants are missing"""
        process = stocal.Process([stocal.Event(['a'], ['b'], 1)])
        traj = self.Sampler(process, {})
        for _ in traj:
            pass
        self.assertEqual(traj.state, stocal.structures.multiset({}))


class TestNextReactionMethod(TestFirstReactionMethod):
    """Test stocal.algorithms.DirectMethod

    This tests the regular TrajectorySampler interface."""
    Sampler = stocal.algorithms.NextReactionMethod


class TestAndersonNRM(TestFirstReactionMethod):
    """Test stocal.algorithms.AndersonNRM"""
    Sampler = stocal.algorithms.AndersonNRM


if __name__ == '__main__':
    unittest.main()
