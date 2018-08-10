"""Bugs reported in github issue tracker
"""
import unittest
import stocal

class Issue1(unittest.TestCase):
    """Issue 1 [CLOSED]: https://github.com/harfel/stocal/issues/1

    Opened 07-02-2018
    Closed 07-02-2018
    """
    def test_text_occurrence_delayed(self):
        """assert that repeated events do not fire before first occurrence"""
        event = stocal.Event(['a'], ['z'], 1, 1)
        self.assertEqual(event.next_occurrence(0), 1)


class Issue2(unittest.TestCase):
    """Issue 2 [CLOSED]: https://github.com/harfel/stocal/issues/2

    Opened 08-02-2018
    Closed 08-02-2018

    traj.transitions is [ 2*a --> aa, 2*b --> bb] with two elements,
    but should be [ 2*a --> aa, 2*a --> aa, a + b --> ab, a + b --> ba,
    2*b --> bb, 2*b --> bb] with six elements.
    """
    def test_pre2017_initial_transitions(self):
        """assert that the initial process has correct potential transitions"""
        from stocal.examples.pre2017 import process
        traj = process.trajectory({'a':10, 'b':10})
        self.assertEqual(len(traj.transitions), 6)


class Issue3(unittest.TestCase):
    """Issue 3 [CLOSED]: https://github.com/harfel/stocal/issues/3

    Opened 21-03-2018
    Closed 21-03-2018
    """
    class Rule(stocal.ReactionRule):
        Transition = stocal.Event

        def novel_reactions(self, x):
            if x=='a':
                yield self.Transition(['a'], [], 10)

    def setUp(self):
        self.process = stocal.Process(transitions=[stocal.Event([], ['a'], 1)],
                                      rules=[self.Rule()])

    def test_expanding_anderson_nrm(self):
        traj = stocal.algorithms.AndersonNRM(self.process, {})
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


class Issue4(unittest.TestCase):
    """Issue 4 [CLOSED]: https://github.com/harfel/stocal/issues/4

    Opened 12-04-2018
    Closed 19-04-2018
    """
    def setUp(self):
        self.process = stocal.Process([stocal.Event(['a'], ['b'], 1)])

    def test_events_with_reactants(self):
        traj = self.process.trajectory({})
        for _ in traj:
            pass
        self.assertEqual(traj.state, stocal.structures.multiset({}))


if __name__ == '__main__':
    unittest.main()
