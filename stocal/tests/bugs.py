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


if __name__ == '__main__':
    unittest.main()
