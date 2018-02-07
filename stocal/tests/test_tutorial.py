"""Test that the exampe is working
"""
import unittest
from stocal import MassAction, Event, ReactionRule, Process


class Dilution(ReactionRule):
    """Dilution rule"""
    Transition = MassAction

    def novel_reactions(self, species):
        yield self.Transition([species], [], 0.001)

class Polymerization(ReactionRule):
    """Polymerization rule"""
    Transition = MassAction

    def novel_reactions(self, k, l):
        yield self.Transition([k, l], [k+l], 10.)

class Hydrolysis(ReactionRule):
    """Hydrolysis rule"""
    Transition = MassAction

    def novel_reactions(self, k):
        for i in range(1, len(k)):
            constant = 10.*i*(len(k)-i)
            yield self.Transition([k], [k[:i], k[i:]], constant)


class TestTutorial(unittest.TestCase):
    """Test every code example of the tutorial"""

    def test_simple_example(self):
        """Starting with a simple system"""
        r1 = MassAction({'A': 2}, {'A2': 1}, 1.)
        r2 = MassAction({'A2': 1}, {'A': 2}, 10.)
        process = Process([r1, r2])
        trajectory = process.trajectory({'A':100}, steps=1000)
        for _ in trajectory:
            result = trajectory.time, trajectory.state.get('A', 0), trajectory.state.get('A2', 0)

    def test_events(self):
        """Adding events"""
        r1 = MassAction({'A': 2}, {'A2': 1}, 1.)
        r2 = MassAction({'A2': 1}, {'A': 2}, 10.)
        feed = Event([], ['A'], 0.0, 1.0)
        process = Process([r1, r2, feed])
        trajectory = process.trajectory({}, steps=100)
        for _ in trajectory:
            pass

    def test_rule_based_dilution(self):
        """Adding dilution"""
        r1 = MassAction({'A': 2}, {'A2': 1}, 1.)
        r2 = MassAction({'A2': 1}, {'A': 2}, 10.)
        feed = Event([], ['A'], 0.0, 1.0)
        process = Process([r1, r2, feed], [Dilution()])
        trajectory = process.trajectory({}, steps=100)
        for _ in trajectory:
            pass

    def test_rule_based_polymers(self):
        """Adding general polymerization and hydrolysis"""
        feed = Event([], ['A'], 0.0, 1.0)
        process = Process(transitions=[feed], rules=[Dilution(), Polymerization(), Hydrolysis()])
        trajectory = process.trajectory({}, steps=100)
        for _ in trajectory:
            pass


if __name__ == '__main__':
    unittest.main()
