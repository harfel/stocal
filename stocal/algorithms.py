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

from .utils import with_metaclass
from .structures import multiset
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
    """
    def __init__(self, process, state, t=0., tmax=float('inf'), steps=None):
        """Initialize the sampler for the given Process and state.

        State is a dictionary that maps chemical species to positive
        integers denoting their copy number. An optional start time
        t, end time tmax, and maximal number of steps can be provided.
        """
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

        for transition in self.process.transitions:
            self.add_transition(transition)
        self.state = multiset()
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
    def propose_transition(self):
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

    def perform_transition(self, time, transition):
        """Perform the given transition.

        Sets sampler.time to time, increases the number of steps,
        performs the given transition by removing reactants and
        adding products to state, and calls Rule.infer_transitions
        for every rule of the process.
        If overwritten by a subclass, the signature can have additional
        arguments which are populated with the argument tuple returned
        by TrajectorySampler.propose_transition.
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

    def has_reached_end(self):
        """True if given max steps or tmax are reached."""
        return self.step == self.steps or self.time >= self.tmax

    def __iter__(self):
        """Sample stochastic trajectory.

        This method iteratively picks an applicable transition by
        calling self.propose_transition and self.perform_transition,
        until one of the given stop criteria are met. During each
        iteration, the firnig transition is yielded.
        """
        while not self.has_reached_end():
            time, transition, args = self.propose_transition()

            if time >= self.tmax:
                break
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
    def __init__(self, process, state, t=0., tmax=float('inf'), steps=None):
        if any(not isinstance(r, Reaction) for r in process.transitions):
            raise ValueError("DirectMethod only works with Reactions.")
        if any(not issubclass(r.Transition, Reaction) for r in process.rules):
            raise ValueError("DirectMethod only works with Reactions.")
        self.transitions = []
        self.propensities = []
        super(DirectMethod, self).__init__(process, state, t, tmax, steps)

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

    def propose_transition(self):
        from random import random
        from math import log

        self.propensities = [r.propensity(self.state) for r in self.transitions]
        total_propensity = sum(self.propensities)
        if not total_propensity:
            return float('inf'), None, tuple()

        delta_t = -log(random())/total_propensity

        transition = None
        pick = random()*total_propensity
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
    def __init__(self, process, state, t=0., tmax=float('inf'), steps=None):
        self.transitions = []
        self.firings = []
        super(FirstReactionMethod, self).__init__(process, state, t, tmax, steps)

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

    def propose_transition(self):
        self.firings = [
            (trans.next_occurrence(self.time, self.state), trans, tuple())
            for trans in self.transitions
        ]

        if self.firings:
            return min(self.firings, key=lambda item: item[0])
        else:
            return float('inf'), None, tuple()


class AndersonNRM(FirstReactionMethod):
    """Next reaction method modified for time-dependent processes

    A stochastic simulation algorithm, published in
    D. F. Anderson, J. Chem. Phys. 127, 214107 (2007)
    as Algorithm (3) 'modified next reaction method'.

    This sampler correctly treats non-autonomous transitions, i.e.
    transitions with time dependent stochastic rates.

    See help(TrajectorySampler) for usage information.
    """
    def __init__(self, process, state, t=0., tmax=float('inf'), steps=None):
        self.T = []
        self.P = []
        super(AndersonNRM, self).__init__(process, state, t, tmax, steps)

    def add_transition(self, transition):
        from math import log
        from random import random

        super(AndersonNRM, self).add_transition(transition)
        self.T.append(0)
        self.P.append(-log(random()))

    def prune_transitions(self):
        depleted = [
            k[0] for t, r, k in self.firings
            if t == float('inf') and (r.rule or isinstance(r, Event))
        ]
        for k in reversed(depleted):
            del self.transitions[k]
            del self.T[k]
            del self.P[k]

    def propose_transition(self):
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
        from random import random

        def int_a_dt(trans, delta_t):
            """Integrate propensity for given delta_t"""
            if isinstance(trans, Event):
                return 0
            else:
                return trans.propensity_integral(self.state, self.time, delta_t)

        self.T = [Tk+int_a_dt(trans, time-self.time) for Tk, trans in
                  zip(self.T, self.transitions)]
        self.P[mu] -= log(random())
        super(AndersonNRM, self).perform_transition(time, transition)
