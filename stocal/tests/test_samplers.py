"""Tests for the stocal.samplers module"""
import unittest
import stocal.samplers


class TestSampler(unittest.TestCase):
    Sampler = stocal.samplers.Sampler

    def test_iter_returns_triple(self):
        """Sampler.__iter__ returns time, state, dict triple"""
        process = stocal.Process([stocal.Event([], ['a'], 1.)])
        sampler = process.sample({})
        result = next(iter(sampler))
        self.assertEqual(len(result), 3)
        self.assertIsInstance(result[0], float)
        self.assertIsInstance(result[1], stocal.multiset)
        self.assertIsInstance(result[2], dict)

    def test_iter_yields_stop_when_empty(self):
        """Sampler.__iter__ raises StopIteration for empty process"""
        process = stocal.Process([])
        with self.assertRaises(StopIteration):
            next(iter(process.sample({})))

    def test_until_time_returns_correct_sampler(self):
        """Sampler.until(time) returns UntilTimeSampler"""
        process = stocal.Process([stocal.Event([], ['a'], 1., 10.)])
        sampler = process.sample({}).until(time=1.)
        self.assertIsInstance(sampler, stocal.samplers.UntilTimeSampler)

    def test_until_steps_returns_correct_sampler(self):
        """Sampler.until(steps) returns UntilStepSampler"""
        process = stocal.Process([stocal.Event([], ['a'], 1., 10.)])
        sampler = process.sample({}).until(steps=10)
        self.assertIsInstance(sampler, stocal.samplers.UntilStepSampler)

    def test_every_time_returns_correct_sampler(self):
        """Sampler.every(time) returns EveryTimeSampler"""
        process = stocal.Process([stocal.Event([], ['a'], 1., 10.)])
        sampler = process.sample({}).every(time=1.)
        self.assertIsInstance(sampler, stocal.samplers.EveryTimeSampler)

    def test_every_steps_returns_correct_sampler(self):
        """Sampler.every(steps) returns EveryStepSampler"""
        process = stocal.Process([stocal.Event([], ['a'], 1., 10.)])
        sampler = process.sample({}).every(steps=10)
        self.assertIsInstance(sampler, stocal.samplers.EveryStepSampler)

    def test_average_time_returns_correct_sampler(self):
        """Sampler.average(time) returns AverageTimeSampler"""
        process = stocal.Process([stocal.Event([], ['a'], 1., 10.)])
        sampler = process.sample({}).average(time=1.)
        self.assertIsInstance(sampler, stocal.samplers.AverageTimeSampler)

    def test_average_steps_returns_correct_sampler(self):
        """Sampler.average(steps) returns AverageStepSampler"""
        process = stocal.Process([stocal.Event([], ['a'], 1., 10.)])
        sampler = process.sample({}).average(steps=10)
        self.assertIsInstance(sampler, stocal.samplers.AverageStepSampler)

    def test_filter_returns_correct_sampler(self):
        """Sampler.filter() returns FilteredSampler"""
        process = stocal.Process([stocal.Event([], ['a'], 1., 10.)])
        sampler = process.sample({}).filter([stocal.MassAction])
        self.assertIsInstance(sampler, stocal.samplers.FilteredSampler)

    def test_state_returns_trajectory_state(self):
        """Sampler.state returns trajectory state"""
        process = stocal.Process([stocal.Event([], ['a'], 1., 10.)])
        state = {}
        sampler = process.sample(state)
        self.assertEqual(sampler.state, state)

    def test_time_returns_trajectory_time(self):
        """Sampler.time returns trajectory time"""
        process = stocal.Process([stocal.Event([], ['a'], 1., 10.)])
        sampler = process.sample({})
        self.assertIsInstance(sampler.time, float)


class TestUntilTimeSampler(TestSampler):
    Sampler = stocal.samplers.UntilTimeSampler

    def test_iter_advances_empty(self, tmax=10.):
        """Sampler.__iter__ advances to end time for empty processes"""
        process = stocal.Process()
        traj = process.sample({})
        sampler = self.Sampler(traj, tmax)
        for _ in sampler:
            pass
        self.assertEqual(sampler.time, tmax)

    def test_iter_advances_long(self, tmax=10.):
        """Sampler.__iter__ advances to end time for long lasting processes"""
        process = stocal.Process([stocal.Event([], ['a'], 1., 100.)])
        traj = process.sample({})
        sampler = self.Sampler(traj, tmax)
        for _ in sampler:
            pass
        self.assertEqual(sampler.time, tmax)

    def test_iter_advances_short(self, tmax=10.):
        """Sampler.__iter__ advances to end time for short lasting processes"""
        process = stocal.Process([stocal.Event(['a'], [], 1., 1.)])
        traj = process.sample({'a': 3})
        sampler = self.Sampler(traj, tmax)
        for _ in sampler:
            pass
        self.assertEqual(sampler.time, tmax)

    def test_iter_includes_all_transitions_at_tmax(self, tmax=1.):
        """Sampler.__iter__ includes all events that happen at tmax"""
        process = stocal.Process([
            stocal.Event([], ['a'], tmax),
            stocal.Event([], ['b'], 0., tmax),
            stocal.Event([], ['c'], tmax/2, tmax/2)])
        traj = process.sample({})
        sampler = self.Sampler(traj, tmax)
        for _ in sampler:
            pass
        self.assertEqual(sampler.step, 5)


class TestUntilStepSampler(TestSampler):
    Sampler = stocal.samplers.UntilStepSampler

    def test_iter_number_of_steps(self, steps=10):
        """Sampler.__iter__ yields exact number of steps"""
        process = stocal.Process([stocal.MassAction([], ['a'], 1.)])
        traj = process.sample({})
        sampler = self.Sampler(traj, steps)
        for _ in sampler:
            pass
        self.assertEqual(sampler.step, steps)

    def test_iter_empty_does_not_proceed(self, steps=10):
        """Sampler.__iter__ does not increase steps for empty process"""
        process = stocal.Process([])
        traj = process.sample({})
        sampler = self.Sampler(traj, steps)
        for _ in sampler:
            pass
        self.assertEqual(sampler.step, 0)


class TestEveryTimeSampler(TestSampler):
    class Sampler(stocal.samplers.EveryTimeSampler):
        def __init__(self, state, skip=False):
            super(TestEveryTimeSampler.Sampler, self).__init__(state, time=1., skip=skip)

    def test_init_optional_skip(self):
        """Sampler.__init__ accepts optional skip argument"""
        sample = stocal.Process().sample({})
        self.Sampler(sample, skip=True)

    def test_iter_yield_times(self):
        """Sampler.__iter__ yields for each interval"""
        process = stocal.Process([stocal.Event([], ['a'], 0.5, 2.4)])
        target = [1., 2., 3., 4., 5.]
        sampler = self.Sampler(process.sample({}))
        for a,b in zip((result[0] for result in sampler), target):
            self.assertEqual(a, b)

    def test_iter_skipping_behavior(self):
        """Sampler.__iter__ skips empty iterations if initialized with skip=True"""
        process = stocal.Process([stocal.Event([], ['a'], 0.5, 2.4)])
        target = [1., 3., 6., 8., 11.]
        sampler = self.Sampler(process.sample({}), skip=True)
        for a,b in zip((result[0] for result in sampler), target):
            self.assertEqual(a, b)

    def test_iter_performs_all_transitions(self, N=20):
        """Sampler.__iter__ performs all transitions"""
        process = stocal.Process([stocal.MassAction(['a'], [], .1)])
        sampler = self.Sampler(process.sample({'a': N}), skip=True)
        total = 0
        for time, state, trans in sampler:
            total += sum(trans.values())
        self.assertEqual(total, N)


class TestEveryStepSampler(TestSampler):
    class Sampler(stocal.samplers.EveryStepSampler):
        def __init__(self, state):
            super(TestEveryStepSampler.Sampler, self).__init__(state, steps=10)

    def test_iter_yield_times(self):
        """Sampler.__iter__ yields for each interval"""
        process = stocal.Process([stocal.Event([], ['a'], 3., 3.)])
        target = [30., 60., 90.]
        sampler = self.Sampler(process.sample({}))
        for a,b in zip((result[0] for result in sampler), target):
            self.assertEqual(a, b)

    def test_iter_performs_all_transitions(self, N=100):
        """Sampler.__iter__ performs all transitions"""
        process = stocal.Process([stocal.MassAction(['a'], [], 1.)])
        sampler = self.Sampler(process.sample({'a': N}))
        total = 0
        for time, state, trans in sampler:
            total += sum(trans.values())
        self.assertEqual(total, N)


class TestAverageTimeSampler(TestEveryTimeSampler):
    class Sampler(stocal.samplers.AverageTimeSampler):
        def __init__(self, state, skip=False):
            super(TestAverageTimeSampler.Sampler, self).__init__(state, time=1., skip=skip)

    def test_iter_correct_averages(self):
        """Sampler.__iter__ calculates correct averages"""
        process = stocal.Process([
            stocal.Event([], ['a'], 0., 3.),
            stocal.Event([], ['a'], 1., 2.)])
        target = [1., 2., 2., 4., 4., 5. ,6., 7., 7., 9., 9., 10., 11.]
        sampler = self.Sampler(process.sample({}))
        for a,b in zip((result[1]['a'] for result in sampler), target):
            self.assertAlmostEqual(a, b)

    def test_iter_correct_skipped_averages(self):
        """Sampler.__iter__ calculates correct averages"""
        process = stocal.Process([
            stocal.Event([], ['a'], 0., 3.),
            stocal.Event([], ['a'], 1., 2.)])
        target = [1., 2., 4., 5. ,6., 7., 9., 10., 11.]
        sampler = self.Sampler(process.sample({}), skip=True)
        for a,b in zip((result[1]['a'] for result in sampler), target):
            self.assertAlmostEqual(a, b)


class TestAverageStepSampler(TestEveryStepSampler):
    class Sampler(stocal.samplers.AverageStepSampler):
        def __init__(self, state, skip=False):
            super(TestAverageStepSampler.Sampler, self).__init__(state, steps=10)

    def test_iter_correct_averages(self):
        """Sampler.__iter__ calculates correct averages"""
        process = stocal.Process([stocal.Event([], ['a'], 0., 3.)])
        target = [5.5, 15.5, 25.5]
        sampler = self.Sampler(process.sample({}))
        for a,b in zip((result[1]['a'] for result in sampler), target):
            self.assertAlmostEqual(a, b)


class TestFilteredSampler(TestSampler):
    class Sampler(stocal.samplers.FilteredSampler):
        def __init__(self, state, transitions=None):
            super(TestFilteredSampler.Sampler, self).__init__(
                state, transitions or [])

    def test_iter_correct_transitions(self):
        """Sampler.__iter__ only yields after declared transitions"""
        r1 = stocal.MassAction(['a'], ['b'], 1.)
        r2 = stocal.MassAction(['b'], ['c'], 1.)
        process = stocal.Process([r1, r2])
        sampler = self.Sampler(process.sample({'a':50, 'b':50}), [r1])
        for result in sampler:
            self.assertIn(r1, result[2])

    def test_iter_empty_filter_list(self):
        """Sampler.__iter__ advances to simulation end when transitions are empty"""
        process = stocal.Process([stocal.MassAction(['a'], [''], 1.)])
        sampler = self.Sampler(process.sample({'a':10}), [])
        with self.assertRaises(StopIteration):
            next(iter(sampler))
        self.assertEqual(sampler.state['a'], 0)


if __name__ == '__main__':
    unittest.main()
