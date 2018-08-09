"""Test that all tutorial examples are working
"""
import unittest
from stocal import MassAction, Event, TransitionRule, Process, multiset
from stocal import algorithms
from stocal.experimental import tauleap

from stocal.tests.test_transitions import TestReactionRule as TestTransitionRule


class Dilution(TransitionRule):
    """Dilution rule"""
    Transition = MassAction

    def novel_reactions(self, species):
        yield self.Transition([species], [], 0.001)

class TestDilution(TestTransitionRule):
    Rule = Dilution


class Polymerization(TransitionRule):
    """Polymerization rule"""
    Transition = MassAction

    def novel_reactions(self, k, l):
        yield self.Transition([k, l], [k+l], 10.)

class TestPolymerization(TestTransitionRule):
    Rule = Polymerization

class Hydrolysis(TransitionRule):
    """Hydrolysis rule"""
    Transition = MassAction

    def novel_reactions(self, k):
        for i in range(1, len(k)):
            constant = 10.*i*(len(k)-i)
            yield self.Transition([k], [k[:i], k[i:]], constant)

class TestHydrolysis(TestTransitionRule):
    Rule = Hydrolysis


class Protein(str):
    pass

class Rna(str):
    pass

class Association(TransitionRule):
    Transition = MassAction
    signature = [Protein, Rna]

    def novel_reactions(self, protein, rna):
        yield self.Transition([protein, rna], [(protein, rna)], 1.)

class TestAssociation(TestTransitionRule):
    Rule = Association


def volume(time, V0=1.0, dV=0.1):
    return V0 + dV*time

class VolumeDependentMassAction(MassAction):
    def propensity(self, state, time):
        a = super(VolumeDependentMassAction, self).propensity(state)
        order = sum(self.reactants.values())
        return a / volume(time)**(order-1)


class TestTutorial(unittest.TestCase):
    """Test every code example of the tutorial"""

    def test_simple_example(self):
        """Starting with a simple system"""
        r1 = MassAction({'A': 2}, {'A2': 1}, 1.)
        r2 = MassAction({'A2': 1}, {'A': 2}, 10.)
        process = Process([r1, r2])
        trajectory = process.sample({'A':100}, steps=1000)
        for _ in trajectory:
            result = trajectory.time, trajectory.state.get('A', 0), trajectory.state.get('A2', 0)

    def test_events(self):
        """Adding events"""
        r1 = MassAction({'A': 2}, {'A2': 1}, 1.)
        r2 = MassAction({'A2': 1}, {'A': 2}, 10.)
        feed = Event([], ['A'], 0.0, 1.0)
        process = Process([r1, r2, feed])
        trajectory = process.sample({}, steps=100)
        for _ in trajectory:
            pass

    def test_rule_based_dilution(self):
        """Adding dilution"""
        r1 = MassAction({'A': 2}, {'A2': 1}, 1.)
        r2 = MassAction({'A2': 1}, {'A': 2}, 10.)
        feed = Event([], ['A'], 0.0, 1.0)
        process = Process([r1, r2, feed], [Dilution()])
        trajectory = process.sample({}, steps=100)
        for _ in trajectory:
            pass

    def test_rule_based_polymers(self):
        """Adding general polymerization and hydrolysis"""
        feed = Event([], ['A'], 0.0, 1.0)
        process = Process(transitions=[feed], rules=[Dilution(), Polymerization(), Hydrolysis()])
        trajectory = process.sample({}, steps=100)
        for _ in trajectory:
            pass

    def test_flattening(self):
        """Flattening a rule-based process"""
        process = Process(rules=[Dilution()])
        flat_process = process.flatten(['a', 'b', 'c'])
        self.assertEqual(len(flat_process.transitions), 3)
        self.assertEqual(len(flat_process.rules), 0)

    def test_types(self):
        """Specifying types via TransitionRule.signature"""
        process = Process(rules=[Association()])
        trajectory = process.sample({Protein('TF'):40,
                                     Rna('mRNA_a'):10,
                                     Rna('mRNA_b'):10},
                                    steps=100)
        flat_process = process.flatten(trajectory.state.domain)
        self.assertEqual(len(list(flat_process.transitions)), 2)        

        for _ in trajectory:
            pass

    def test_time_dependence(self):
        """Specifying volume-dependent reactions"""
        process = Process([VolumeDependentMassAction(['x', 'x'], ['x2'], 1.)])
        trajectory = process.sample({'x':100}, steps=100)
        for _ in trajectory:
            pass

    def test_samplers_available(self):
        self.assertTrue(issubclass(algorithms.DirectMethod,
                                   algorithms.StochasticSimulationAlgorithm))
        self.assertTrue(issubclass(algorithms.FirstReactionMethod,
                                   algorithms.StochasticSimulationAlgorithm))
        self.assertTrue(issubclass(algorithms.NextReactionMethod,
                                   algorithms.StochasticSimulationAlgorithm))
        self.assertTrue(issubclass(algorithms.AndersonMethod,
                                   algorithms.StochasticSimulationAlgorithm))
        self.assertTrue(issubclass(tauleap.CaoMethod,
                                   algorithms.StochasticSimulationAlgorithm))


if __name__ == '__main__':
    unittest.main()
