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
    imlementations of the stochastic simulation algorithm.
    A trajectory sampler is initialized with a given process and
    initial state, and options start time, end time, and maximal
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

        self.transitions = []
        for transition in self.process.transitions:
            self.add_transition(transition)
        self.state = multiset()
        self.update_state(state)

    @abc.abstractmethod
    def add_transition(self, transition):
        """Add a new transition to the sampler"""
        pass

    @abc.abstractmethod
    def update_state(self, dct, **kwds):
        """Modify sampler state.

        New transitions get infered appropriately."""
        pass

    @abc.abstractmethod
    def __iter__(self):
        """Sample stochastic trajectory.

        This method iteratively picks an applicable transition with
        probability proportional to its propensity, updates the sampler
        state accordingly, and yields the transition to the caller.
        After each transition, rules are invoked to infer potential
        novel transitions.
        """
        raise StopIteration

    def _perform_transition(self, transition):
        """Perform the given transition.

        Remove reactants from and add products to state, and
        call self.rules to potentially infer novel transitions
        for the changed system state.
        """
        self.state -= transition.true_reactants
        for rule in self.process.rules:
            # XXX determine from state which rules should be queried
            for trans in rule.infer_transitions(transition.true_products, self.state):
                trans.rule = rule
                self.add_transition(trans)
        self.state += transition.true_products


class DirectMethod(TrajectorySampler):
    """Implementation of Gillespie's direct method.

    The standard stochastic simulation algorithm, published in
    D. T. Gillespie, J. Comp. Phys. 22, 403-434 (1976).

    DirectMethod works only over processes that employ stochastic
    transitions (Reaction's). If the process includes Event's,
    FirstReactionMethod has to be used instead.

    See help(TrajectorySampler) for usage information.
    """
    def __init__(self, process, state, t=0., tmax=float('inf'), steps=None):
        super(DirectMethod, self).__init__(process, state, t, tmax, steps)
        if any(not isinstance(r, Reaction) for r in process.transitions):
            raise ValueError("DirectMethod only works with Reactions.")
        if any(not issubclass(r.Transition, Reaction) for r in process.rules):
            raise ValueError("DirectMethod only works with Reactions.")

    def add_transition(self, transition):
        self.transitions.append(transition)

    def update_state(self, dct):
        for rule in self.process.rules:
            for trans in rule.infer_transitions(dct, self.state):
                if trans not in self.transitions:
                    trans.rule = rule
                    self.add_transition(trans)
        self.state.update(dct)

    def __iter__(self):
        from random import random
        from math import log

        while True:
            if self.steps is not None and self.step == self.steps:
                return
            if self.time >= self.tmax:
                break

            propensities = [r.propensity(self.state) for r in self.transitions]

            total_propensity = sum(propensities)
            if not total_propensity:
                break

            # housekeeping: remove depleted transitions
            depleted = [
                i for i, (p, t) in enumerate(zip(propensities, self.transitions))
                if p == 0. and t.rule
            ]
            for i in reversed(depleted):
                del self.transitions[i]
                del propensities[i]

            delta_t = -log(random())/total_propensity
            self.time += delta_t

            if self.time >= self.tmax:
                break
            else:
                transition = None
                pick = random()*total_propensity
                for propensity, transition in zip(propensities, self.transitions):
                    pick -= propensity
                    if pick < 0.:
                        break
                self.step += 1
                self._perform_transition(transition)
                yield transition

        if self.tmax < float('inf'):
            self.time = self.tmax


class FirstReactionMethod(TrajectorySampler):
    """Implementation of Gillespie's first reaction method.

    A stochastic simulation algorithm, published in
    D. T. Gillespie, J. Comp. Phys. 22, 403-434 (1976).

    FirstReactionMethod works with processes that feature deterministic
    transitions, i.e. Event's.

    See help(TrajectorySampler) for usage information.
    """
    def __iter__(self):
        while True:
            if self.steps is not None and self.step == self.steps:
                return
            if self.time >= self.tmax:
                break

            firings = [
                (trans.next_occurrence(self.time, self.state), trans)
                for trans in self.transitions
            ]

            if not firings:
                break
            time, transition = min(firings, key=lambda item: item[0])

            # housekeeping: remove depleted transitions
            depleted = [
                i for i, (t, r) in enumerate(firings)
                if t == float('inf') and (r.rule or isinstance(r, Event))
            ]
            for i in reversed(depleted):
                del self.transitions[i]
                del firings[i]

            if time >= self.tmax:
                break
            else:
                self.step += 1
                self.time = time
                transition.last_occurrence = time
                self._perform_transition(transition)
                yield transition

        if self.tmax < float('inf'):
            self.time = self.tmax

    def add_transition(self, transition):
        self.transitions.append(transition)

    def update_state(self, dct):
        for rule in self.process.rules:
            for trans in rule.infer_transitions(dct, self.state):
                if trans not in self.transitions:
                    trans.rule = rule
                    self.add_transition(trans)
        self.state.update(dct)
