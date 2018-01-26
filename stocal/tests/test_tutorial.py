import unittest
from stocal import MassAction, Event, ReactionRule, Process


class Dilution(ReactionRule) :
    Transition = MassAction
    order = 1

    def novel_reactions(self, species) :
        yield self.Transition([species], [], 0.001)

class Polymerization(ReactionRule) :
    Transition = MassAction
    order = 2

    def novel_reactions(self, k, l) :
        yield self.Transition([k,l], [k+l], 10.)

class Hydrolysis(ReactionRule) :
    Transition = MassAction
    order = 1

    def novel_reactions(self, k) :
        for i in xrange(1, len(k)) :
            c = 10.*i*(len(k)-i)
            yield self.Transition([k], [k[:i], k[i:]], c)


class TestTutorial(unittest.TestCase) :
	def test_simple_example(self) :
		r1 = MassAction({'A': 2}, {'A2': 1}, 1.)
		r2 = MassAction({'A2': 1}, {'A': 2}, 10.)
		process = Process([r1, r2])
		trajectory = process.trajectory({'A':100}, steps=1000)
		for transition in trajectory :
			trajectory.time, trajectory.state.get('A', 0), trajectory.state.get('A2', 0)

	def test_events(self) :
		r1 = MassAction({'A': 2}, {'A2': 1}, 1.)
		r2 = MassAction({'A2': 1}, {'A': 2}, 10.)
		feed = Event([], ['A'], 0.0, 1.0)
		process = Process([r1, r2, feed])
		trajectory = process.trajectory({}, steps=100)
		for transition in trajectory : pass

	def test_rule_based_dilution(self) :
		r1 = MassAction({'A': 2}, {'A2': 1}, 1.)
		r2 = MassAction({'A2': 1}, {'A': 2}, 10.)
		feed = Event([], ['A'], 0.0, 1.0)
		process = Process([r1, r2, feed], [Dilution()])
		trajectory = process.trajectory({}, steps=100)
		for transition in trajectory : pass

	def test_rule_based_polymers(self) :
		feed = Event([], ['A'], 0.0, 1.0)
		process = Process(transitions=[feed], rules=[Dilution(), Polymerization(), Hydrolysis()])
		trajectory = process.trajectory({}, steps=100)
		for transition in trajectory : pass
