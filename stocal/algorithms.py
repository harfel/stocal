# Copyright 2018 Harold Fellermann
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""Stochastic simulation algorithms

This module provides implementations of different stochastic simulation
algorithms. All implementations have to implement the TrajectorySampler
interface by subclassing this abstract base class.
"""

import abc

from ._utils import with_metaclass
from .structures import multiset, DependencyGraph, QueueWrapper
from .transitions import Event, Reaction


class TrajectorySampler(with_metaclass(abc.ABCMeta, object)):
    """Abstract base class for stochastic trajectory samplers.

    This is the interface for stochastic trajectory samplers, i.e.
    implementations of the stochastic simulation algorithm.
    A trajectory sampler is initialized with a given process and
    initial state, and optional start time, end time, and maximal
    number of iterations. The sampler instance can then be iterated
    over to produce a stochastic trajectory:

    >>> trajectory = DirectMethod(process, state, steps=10000)
    >>> for transition in trajectory:
    ...     print trajectory.time, trajectory.state, transition

    When implementing a novel TrajectorySampler, make sure to only
    generate random numbers using the TrajectorySampler.rng random
    number generator.
    """
    def __init__(self, process, state, t=0., tmax=float('inf'), steps=None, seed=None):
        """Initialize the sampler for the given Process and state.

        State is a dictionary that maps chemical species to positive
        integers denoting their copy number. An optional start time
        t, end time tmax, maximal number of steps, and random seed
        can be provided.
        """
        from random import Random

        if t < 0:
            raise ValueError("t must not be negative.")
        if tmax < 0:
            raise ValueError("tmax must not be negative.")
        if steps is not None and steps < 0:
            raise ValueError("steps must not be negative.")
        if any(n <= 0 for n in state.values()):
            raise ValueError("state copy numbers must be positive")

        self.process = process
        self.step = 0
        self.steps = steps
        self.time = t
        self.tmax = tmax
        self.rng = Random(seed)

        self.transitions = []
        self.state = multiset()
        self.update_state(state)

        for transition in self.process.transitions:
            self.add_transition(transition)

        self.update_state(state)

    def update_state(self, dct):
        """Modify sampler state.

        Update system state and infer new applicable transitions.
        """
        for rule in self.process.rules:
            for trans in rule.infer_transitions(dct, self.state):
                trans.rule = rule
                self.add_transition(trans)
        self.state.update(dct)

    @abc.abstractmethod
    def add_transition(self, transition):
        """Add a new transition to the sampler

        Must be implemented by a subclass.
        """
        pass

    @abc.abstractmethod
    def prune_transitions(self):
        """Remove infered transitions that are no longer applicable.

        Must be implemented by a subclass.
        """
        pass

    @abc.abstractmethod
    def propose_potential_transition(self):
        """Propose new transition

        Must be implemented by a subclass.
        Must return a triple where the first item is the time of the
        proposed transition, the second item the transition itself,
        and the last item a tuple of additional arguments that are
        passed through to calls to TrajectorySampler.perform_transition.
        Time and transition have to be picked from the correct
        probability distributions.
        """
        return float('inf'), None, tuple()

    def is_applicable(self, time, transition, *args):
        """True if the transition is applicable
        
        The standard implementation always returns True, i.e. it assumes
        that any transition returned by propose_potential_transition is
        applicable. Overwrite this method if you want to implement
        accept/reject type of samplers such as composition-rejection.
        """
        return True

    def perform_transition(self, time, transition):
        """Perform the given transition.

        Sets sampler.time to time, increases the number of steps,
        performs the given transition by removing reactants and
        adding products to state, and calls Rule.infer_transitions
        for every rule of the process.
        If overwritten by a subclass, the signature can have additional
        arguments which are populated with the argument tuple returned
        by TrajectorySampler.propose_potential_transition.
        """
        self.step += 1
        self.time = time
        transition.last_occurrence = time
        self.state -= transition.true_reactants
        for rule in self.process.rules:
            for trans in rule.infer_transitions(transition.true_products, self.state):
                trans.rule = rule
                self.add_transition(trans)
        self.state += transition.true_products
        self.prune_transitions()

    def reject_transition(self, time, transition, *args):
        """Do not execute the given transition.
        
        The default implementation does nothing. Overwrite this method
        if you, for example, want to prevent the same transition from
        being proposed again.
        """
        pass

    def has_reached_end(self):
        """True if given max steps or tmax are reached."""
        return self.step == self.steps or self.time >= self.tmax

    def __iter__(self):
        """Standard interface to sample a stochastic trajectory.

        Yields each performed transition of the stochastic trajectory.

        This implementation first picks a potential transition.
        If the transition is applicable, it is performed, otherwise
        it is rejected. Iteration continues until a stop criterion
        occurs.
        Consider to overwrite self.propose_potential_transition,
        self.is_applicable, self.perform_transition,
        self.reject_transition or self.has_reached_end in favour of
        overloading __iter__ when implementing a TrajectorySampler.
        """
        while not self.has_reached_end():
            time, transition, args = self.propose_potential_transition()

            if time >= self.tmax:
                break
            elif not self.is_applicable(time, transition, *args):
                self.reject_transition(time, transition, *args)
            else:
                self.perform_transition(time, transition, *args)
                yield transition

        if self.step != self.steps and self.tmax < float('inf'):
            self.time = self.tmax


class DirectMethod(TrajectorySampler):
    """Implementation of Gillespie's direct method.

    The standard stochastic simulation algorithm, published in
    D. T. Gillespie, J. Comp. Phys. 22, 403-434 (1976).

    DirectMethod works only over processes that employ stochastic
    transitions whose propensity functions do not explicitly depend on
    time (autonomous Reaction's).

    DirectMethod maintains a list of current transitions and a list
    of propensities. When proposing transitions, propensities are
    calculated for all transitions and time and occuring transition
    are drawn from the appropriate probability distributions.

    See help(TrajectorySampler) for usage information.
    """
    def __init__(self, process, state, t=0., tmax=float('inf'), steps=None, seed=None):
        if any(not isinstance(r, Reaction) for r in process.transitions):
            raise ValueError("DirectMethod only works with Reactions.")
        if any(not issubclass(r.Transition, Reaction) for r in process.rules):
            raise ValueError("DirectMethod only works with Reactions.")
        self.transitions = []
        self.propensities = []
        super(DirectMethod, self).__init__(process, state, t, tmax, steps, seed)

    def add_transition(self, transition):
        self.transitions.append(transition)

    def prune_transitions(self):
        depleted = [
            i for i, (p, t) in enumerate(zip(self.propensities, self.transitions))
            if p == 0. and t.rule
        ]
        for i in reversed(depleted):
            del self.transitions[i]
            del self.propensities[i]

    def propose_potential_transition(self):
        from math import log

        self.propensities = [r.propensity(self.state) for r in self.transitions]
        total_propensity = sum(self.propensities)
        if not total_propensity:
            return float('inf'), None, tuple()

        delta_t = -log(self.rng.random())/total_propensity

        transition = None
        pick = self.rng.random()*total_propensity
        for propensity, transition in zip(self.propensities, self.transitions):
            pick -= propensity
            if pick < 0.:
                break

        return self.time + delta_t, transition, tuple()


class FirstReactionMethod(TrajectorySampler):
    """Implementation of Gillespie's first reaction method.

    A stochastic simulation algorithm, published in
    D. T. Gillespie, J. Comp. Phys. 22, 403-434 (1976).

    FirstReactionMethod works with processes that feature deterministic
    transitions, i.e. Event's.

    See help(TrajectorySampler) for usage information.
    """
    def __init__(self, process, state, t=0., tmax=float('inf'), steps=None, seed=None):
        self.transitions = []
        self.firings = []
        super(FirstReactionMethod, self).__init__(process, state, t, tmax, steps, seed)

    def add_transition(self, transition):
        self.transitions.append(transition)

    def prune_transitions(self):
        depleted = [
            i for i, (t, r, _) in enumerate(self.firings)
            if t == float('inf') and (r.rule or isinstance(r, Event))
        ]
        for i in reversed(depleted):
            del self.transitions[i]
            del self.firings[i]

    def propose_potential_transition(self):
        self.firings = [
            (trans.next_occurrence(self.time, self.state, self.rng), trans, tuple())
            for trans in self.transitions
        ]

        if self.firings:
            return min(self.firings, key=lambda item: item[0])
        else:
            return float('inf'), None, tuple()

    def is_applicable(self, time, transition, *args):
        """Returns False for Event's that lack their reactants."""
        if isinstance(transition, Event):
            return transition.reactants <= self.state
        else :
            return super(FirstReactionMethod, self).is_applicable(time, transition, *args)

    def reject_transition(self, time, transition, *args):
        """Reject inapplicable Event
        
        Advance system time to transition time and pretend transition
        had happened there, but do not change the state.
        """
        self.time = time
        transition.last_occurrence = time


class AndersonNRM(FirstReactionMethod):
    """Next reaction method modified for time-dependent processes

    A stochastic simulation algorithm, published in
    D. F. Anderson, J. Chem. Phys. 127, 214107 (2007)
    as Algorithm (3) 'modified next reaction method'.

    This sampler correctly treats non-autonomous transitions, i.e.
    transitions with time dependent stochastic rates.

    See help(TrajectorySampler) for usage information.
    """
    def __init__(self, process, state, t=0., tmax=float('inf'), steps=None, seed=None):
        self.T = []
        self.P = []
        super(AndersonNRM, self).__init__(process, state, t, tmax, steps, seed)

    def add_transition(self, transition):
        from math import log

        super(AndersonNRM, self).add_transition(transition)
        self.T.append(0)
        self.P.append(-log(self.rng.random()))

    def prune_transitions(self):
        depleted = [
            k[0] for t, r, k in self.firings
            if t == float('inf') and (r.rule or isinstance(r, Event))
        ]
        for k in reversed(depleted):
            del self.transitions[k]
            del self.T[k]
            del self.P[k]

    def propose_potential_transition(self):
        def eq_13(trans, target):
            """Determine timestep for a transition in global time scale"""
            if isinstance(trans, Event):
                return trans.next_occurrence(self.time)
            else:
                return trans.propensity_meets_target(self.state, self.time, target) +  self.time

        self.firings = [
            (eq_13(trans, Pk-Tk), trans, (k,))
            for k, (trans, Pk, Tk)
            in enumerate(zip(self.transitions, self.P, self.T))
        ]
        if self.firings:
            return min(self.firings, key=lambda item: item[0])
        else:
            return float('inf'), None, tuple()

    def perform_transition(self, time, transition, mu):
        from math import log

        def int_a_dt(trans, delta_t):
            """Integrate propensity for given delta_t"""
            if isinstance(trans, Event):
                return 0
            else:
                return trans.propensity_integral(self.state, self.time, delta_t)

        self.T = [Tk+int_a_dt(trans, time-self.time) for Tk, trans in
                  zip(self.T, self.transitions)]
        self.P[mu] -= log(self.rng.random())
        super(AndersonNRM, self).perform_transition(time, transition)


class NextReactionMethod(TrajectorySampler):
    def __init__(self, process, state, t=0., tmax=float('inf'), steps=None, seed=None):
        self.dependency_graph = DependencyGraph(process.transitions)
        self.queue_wrapper = QueueWrapper()
        super(NextReactionMethod, self).__init__(process, state, t, tmax, steps, seed)
        self.queue_wrapper.initialise_transitions(self.time, self.state, self.rng)

    def propose_potential_transition(self):
        if self.queue_wrapper.queue:
            trans, time = self.queue_wrapper.queue.topitem()
            return time, trans[0], trans[1:]
        else:
            return float('inf'), None, tuple()

    def is_applicable(self, time, transition, *args):
        if isinstance(transition, Event):
            return transition.reactants <= self.state
        else :
            return super(NextReactionMethod, self).is_applicable(time, transition, *args)

    def perform_transition(self, time, transition, mult=0):
        super(NextReactionMethod, self).perform_transition(time, transition)
        self.queue_wrapper.update_transitions(transition, mult, time, self.state, self.dependency_graph, self.rng)

    def reject_transition(self, time, transition, mult):
        self.time = time
        transition.last_occurrence = time
        self.queue_wrapper.update_transitions(transition, mult, time, self.state, self.dependency_graph, self.rng)

    def add_transition(self, transition):
        self.transitions.append(transition)
        self.dependency_graph.add_reaction(transition)
        self.queue_wrapper.add_transition(transition, self.time, self.state)

    def remove_transition(self, transition, multiplicity):
        self.transitions.remove(transition)
        self.dependency_graph.remove_reaction(transition)
        self.queue_wrapper.remove_transition(transition, multiplicity)

    def update_state(self, dct):
        for rule in self.process.rules:
            for trans in rule.infer_transitions(dct, self.state):
                trans.rule = rule
                self.add_transition(trans)
        self.state.update(dct)
        self.queue_wrapper.update_state_transitions(dct, self.time, self.state, self.dependency_graph, self.rng)

    def prune_transitions(self):
        for trans, time in self.queue_wrapper.queue.items():
            if time == float('inf') and (trans[0].rule or isinstance(trans[0], Event)):
                self.remove_transition(*trans)
