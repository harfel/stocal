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
algorithms. All implementations have to implement the
StochasticSimulationAlgorithm interface by subclassing this abstract
base class.

The module also provides several helper classes for general stochastic
simulation algorithms. DependencyGraph is a generic species-dependency
graph that gives quick access to the set of affected transitions from
a set of affected species. MultiDict and PriorityQueue are data
structures that can be used for managing transitions. MultiDict is
a mapping from transitions to propensities used in DirectMethod and
related algorithms. Priority Queue is an indexed priority queue as it
appears in Gibson-Bruck like NextReactionMethod's. Both data structures
allow for transitions to be added multiple times, and keep track of
their multiplicity.
"""

import sys
import abc
import logging
import warnings
from pqdict import pqdict

from ._utils import with_metaclass
from .structures import multiset
from .transitions import Event, Reaction


class DependencyGraph(dict):
    """Species-transition dependency graph

    A mapping from species to transitions that are affected by a
    change in the respecive species' count. A transition counts as
    affected if it appears among the reactants of the transition.
    """
    def add_reaction(self, reaction):
        """Add reaction to dependencies"""
        for reactant in reaction.reactants:
            self.setdefault(reactant, set()).add(reaction)

    def remove_reaction(self, reaction):
        """Remove reaction from dependencies"""
        for reactant in reaction.affected_species:
            if reactant in self:
                self[reactant].discard(reaction)
                if not self[reactant]:
                    del self[reactant]

    def affected_transitions(self, species):
        """Get all transitions affected by a change in any of the given species

        Species is an iterable and the method returns a set of Transition's.
        """
        return set(
            trans for reactant in species
            for trans in self.get(reactant, [])
        )


class MultiDict(object):
    """Dictionary with multiplicity count

    MultiDict supports DirectMethod like samplers with a dictionary
    that maps transitions to propensities. Unlike normal dictionaries
    MultiDict keep count of the number of times a key is added, i.e.
    its multipliciry. Therefore, the signature of MultiDict is

        transition -> (propensity, multiplicity)

    This class and its interface are not module level implementation
    details and are not specified through tests. They might change in
    future versions of stocal.
    """
    def __init__(self):
        self._dict = dict()

    def __contains__(self, item):
        return item not in self._dict

    def add_item(self, key, value):
        """Add item with associated value.

        If the key has been inserted ealrier, its multiplicity is
        increased by one. Otherwise, it is stored with multiplicity
        one."""
        if key not in self._dict:
            # insert transition with propensity and multiplicity one
            self._dict[key] = [value, 1]
        else:
            # increase multiplicity count
            self._dict[key][1] += 1

    def __delitem__(self, key):
        """Delete all instances of the given key."""
        del self._dict[key]

    def keys(self):
        """Return a list of all present transition instances."""
        return (
            key for key, (prop, mult) in self._dict.items()
            for i in range(mult)
        )

    def items(self):
        """Iterate over key, value, multiplicity triples."""
        for key, item in self._dict.items():
            value, multiplicity = item
            yield key, value, multiplicity

    def update_item(self, key, value):
        """Change value associated with key"""
        self._dict[key][0] = value


class PriorityQueue(object):
    """Indexed priority queue

    Data structure used for Gibson-Bruck-like transition selection.
    See the documentation of pqdict for the general properties of the
    data structure. Unlike the standard indexed priority queue, this
    implementation allows keys (transitions) to have multiple associated
    values. In addition, each instance of a key can have optional
    associated data. The exact mapping is therefore

        key -> [(value_1, data_1), (value_2, data_2), ...]


    Instead of assigning values directly, The PriorityQueue constructor
    takes a function to calculate values for a key. Its signature is

        value_callback(key, data) -> float.

    Custom data is passed as keyword arguments to the add_item method.
    add_item creates a mutable instance of type PriorityQueue.Item for
    each transition and binds keyword arguments to it. This data object
    is retrieved by queue access methods and is also passed to the
    value_callback.
    """
    class Item(object):
        """Namespace to hold custom associated data"""
        def __init__(self, **params):
            for key, value in params.items():
                setattr(self, key, value)

        def __repr__(self):
            return '<Data %s>' % ', '.join('%s=%r' % (k, v)
                                           for k, v
                                           in self.__dict__.items())

        def __eq__(self, other):
            return self.__dict__ == other.__dict__

        def __lt__(self, other):
            return False


    def __init__(self, value_callback):
        self._queue = pqdict()
        self.value_callback = value_callback
        self.depletion_value = float('inf')

    def __bool__(self):
        return bool(self._queue)

    __nonzero__ = __bool__

    def __getitem__(self, key):
        return self._queue[key]

    def keys(self):
        """Return a list of all present transition instances."""
        return (
            key for key, instances in self._queue.items()
            for i in instances
        )

    def topitem(self):
        """Yield (value, key, (data,)) triple with lowest value."""
        key, instances = self._queue.topitem()
        value, data = instances[0]
        return value, key, (data,)

    def add_item(self, key, **params):
        """Add instance of key to the queue.

        Any additional keyword arguments are used to populate a
        PriorityQueue.Item object that is subsequently passed to
        the queue's value_callback and stored together with the value
        of the instance.
        """
        if key not in self._queue:
            self._queue[key] = []

        data = self.Item(**params)
        value = self.value_callback(key, data)
        self._queue[key].append((value, data))
        self._queue[key].sort()
        self._queue.heapify()

    def remove_item(self, key):
        """Remove all instances of key with a value of float('inf')"""
        # XXX The code below is wrong and leads to a memory leak. The section should read
        # remaining = [t for t in self._queue[key]
        #             if t[0] != self.depletion_value]
        # but this leads to an inconsistent sampler state (orphaned reactions in the
        # dependency graph). The fault needs to be fixed, but the error is not only at
        # this point.
        remaining = [t for t in self._queue[key]
                     if t != self.depletion_value]
        if remaining:
            self._queue[key] = remaining
        else:
            del self._queue[key]

    def update_one_instance(self, key):
        """Recalculates the value for the 'first' instance of key."""
        instances = self._queue[key]
        _, data = instances[0]
        instances[0] = self.value_callback(key, data), data
        instances.sort()
        self._queue.heapify()

    def update_items(self, keys):
        """Recalculate values for all provided keys."""
        for key in keys:
            instances = sorted(
                (self.value_callback(key, data), data)
                for value, data in self._queue[key]
            )
            self._queue[key] = instances


class StochasticSimulationAlgorithm(with_metaclass(abc.ABCMeta, object)):
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

    When implementing a novel StochasticSimulationAlgorithm, make sure
    to only generate random numbers using the
    StochasticSimulationAlgorithm.rng random number generator.
    """
    def __init__(self, process, state, tstart=0., tmax=float('inf'), steps=None, seed=None):
        """Initialize the sampler for the given Process and state.

        State is a dictionary that maps chemical species to positive
        integers denoting their copy number. An optional start time
        t, end time tmax, maximal number of steps, and random seed
        can be provided.
        """
        from random import Random

        if tstart < 0:
            raise ValueError("t must not be negative.")
        if tmax < 0:
            raise ValueError("tmax must not be negative.")
        if steps is not None and steps < 0:
            raise ValueError("steps must not be negative.")
        if any(n < 0 for n in state.values()):
            raise ValueError("state copy numbers must be positive")

        self.process = process
        self.step = 0
        self.steps = steps
        self.time = tstart
        self.tmax = tmax
        self.rng = Random(seed)

        self.state = multiset()
        self.update_state(state)

        for transition in self.process.transitions:
            self.add_transition(transition)

    @abc.abstractproperty
    def transitions(self):
        """Return list of all transitions

        Implementations need to return every instance of an added
        transition, i.e. two copies of a transition that has
        multiplicity two.
        """
        return []

    def update_state(self, mapping):
        """Modify sampler state.

        Update system state and infer new applicable transitions.
        """
        mapping = mapping if isinstance(mapping, multiset) else multiset(mapping)
        for rule in self.process.rules:
            for trans in rule.infer_transitions(mapping, self.state):
                trans.rule = rule
                self.add_transition(trans)
        self.state.update(mapping)

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
        passed through to calls to
        StochasticSimulationAlgorithm.perform_transition.
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

    def perform_transition(self, time, transition, *args):
        """Perform the given transition.

        Sets sampler.time to time, increases the number of steps,
        performs the given transition by removing reactants and
        adding products to state, and calls Rule.infer_transitions
        for every rule of the process.
        If overwritten by a subclass, the signature can have additional
        arguments which are populated with the argument tuple returned
        by StochasticSimulationAlgorithm.propose_potential_transition.
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

    def reject_transition(self, time, transition, *args):
        """Do not execute the given transition.

        The default implementation does nothing. Overwrite this method
        if you, for example, want to prevent the same transition from
        being proposed again.
        """
        pass

    def has_reached_end(self):
        """True if given max steps or tmax are reached."""
        return (self.steps is not None and self.step >= self.steps) or self.time > self.tmax

    def __iter__(self):
        return self

    def __next__(self):
        """Interface to sample a stochastic trajectory.

        Performs a single random step of the Gillespie algorithm and
        returns the fired transition after updating self.time and
        self.state and inferring potentially novel transitions.
        Consider to overwrite self.propose_potential_transition,
        self.is_applicable, self.perform_transition,
        self.reject_transition or self.has_reached_end in favour of
        overloading __next__ when implementing a
        StochasticSimulationAlgorithm.
        """
        # XXX can this logic be simplified?
        while not self.has_reached_end():
            time, transition, args = self.propose_potential_transition()

            if time > self.tmax or time == float('inf'):
                break
            elif not self.is_applicable(time, transition, *args):
                self.reject_transition(time, transition, *args)
            else:
                self.perform_transition(time, transition, *args)
                self.prune_transitions()
                return transition

        if self.step != self.steps and self.tmax < float('inf'):
            self.time = self.tmax
        raise StopIteration

    def until(self, time):
        old_tmax, self.tmax = self.tmax, time
        while self.time < time:
            try:
                next(self)
            except StopIteration:
                break
        self.tmax = old_tmax
        # XXX what to return?

class DirectMethod(StochasticSimulationAlgorithm):
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

    See help(StochasticSimulationAlgorithm) for usage information.
    """
    def __init__(self, process, state, **sampler_options):
        if any(not isinstance(r, Reaction) for r in process.transitions):
            raise ValueError("DirectMethod only works with Reactions.")
        if any(not issubclass(r.Transition, Reaction) for r in process.rules):
            raise ValueError("DirectMethod only works with Reactions.")
        self.dependency_graph = DependencyGraph()
        self.propensities = MultiDict()
        self.depleted = []
        super(DirectMethod, self).__init__(process, state, **sampler_options)

    @property
    def transitions(self):
        return self.propensities.keys()

    def update_state(self, dct):
        super(DirectMethod, self).update_state(dct)
        affected_transitions = self.dependency_graph.affected_transitions(dct)
        self.update_propensities(affected_transitions)

    def add_transition(self, transition):
        self.dependency_graph.add_reaction(transition)
        propensity = transition.propensity(self.state)
        self.propensities.add_item(transition, propensity)

    def prune_transitions(self):
        for trans in self.depleted:
            self.dependency_graph.remove_reaction(trans)
            del self.propensities[trans]
        self.depleted = []

    def propose_potential_transition(self):
        from math import log

        total_propensity = sum(mult*prop
                               for transition, prop, mult
                               in self.propensities.items())
        if not total_propensity:
            return float('inf'), None, tuple()

        delta_t = -log(self.rng.random())/total_propensity

        transition = None
        pick = self.rng.random()*total_propensity
        for transition, prop, mult in self.propensities.items():
            pick -= mult*prop
            if pick < 0.:
                break

        return self.time + delta_t, transition, tuple()

    def perform_transition(self, time, transition):
        super(DirectMethod, self).perform_transition(time, transition)
        affected = self.dependency_graph.affected_transitions(transition.affected_species)
        self.update_propensities(affected)

    def update_propensities(self, affected_transitions):
        """Update propensities of all given transitions.

        If the new propensity of a transition is 0 and the transition
        has been derived by a rule or is an Event, the transition gets
        added to self.depleted for later pruning."""
        # XXX order in which reactions are updated needs to be preserved over independent runs
        for trans in affected_transitions:
            propensity = trans.propensity(self.state)
            if propensity == 0 and (trans.rule or isinstance(trans, Event)):
                self.depleted.append(trans)
            self.propensities.update_item(trans, propensity)


class FirstReactionMethod(StochasticSimulationAlgorithm):
    """Implementation of Gillespie's first reaction method.

    A stochastic simulation algorithm, published in
    D. T. Gillespie, J. Comp. Phys. 22, 403-434 (1976).

    FirstReactionMethod works with processes that feature deterministic
    transitions, i.e. Event's.

    See help(StochasticSimulationAlgorithm) for usage information.
    """
    def __init__(self, process, state, **sampler_options):
        self._transitions = []
        self.firings = []
        super(FirstReactionMethod, self).__init__(process, state, **sampler_options)

    @property
    def transitions(self):
        return self._transitions

    def add_transition(self, transition):
        self._transitions.append(transition)

    def prune_transitions(self):
        depleted = [
            i for i, (t, r, _) in enumerate(self.firings)
            if t == float('inf') and (r.rule or isinstance(r, Event))
        ]
        for i in reversed(depleted):
            del self._transitions[i]
            del self.firings[i]

    def propose_potential_transition(self):
        self.firings = [
            (trans.next_occurrence(self.time, self.state, self.rng), trans, tuple())
            for trans in self._transitions
        ]

        if self.firings:
            return min(self.firings, key=lambda item: item[0])
        else:
            return float('inf'), None, tuple()

    def is_applicable(self, time, transition, *args):
        """Returns False for Event's that lack their reactants."""
        if isinstance(transition, Event):
            return transition.reactants <= self.state
        else:
            return super(FirstReactionMethod, self).is_applicable(time, transition, *args)

    def reject_transition(self, time, transition, *args):
        """Reject inapplicable Event

        Advance system time to transition time and pretend transition
        had happened there, but do not change the state.
        """
        self.time = time
        transition.last_occurrence = time


class NextReactionMethod(FirstReactionMethod):
    """Implementation of Gibson & Bruck's next reaction method.

    A stochastic simulation algorithm, published in
    M. A. Gibson & J. Bruck, J. Phys. Chem. A 2000, 104, 1876-1889 (1999).

    NextReactionMethod works with processes that feature deterministic
    transitions, i.e. Event's. It is an optimized version of Gillespie's
    FirstReactionMethod whose runtime scales logarithmically with the
    number of reactions.

    See help(StochasticSimulationAlgorithm) for usage information.
    """
    def __init__(self, process, state, **sampler_options):
        self.dependency_graph = DependencyGraph()
        self.firings = PriorityQueue(self.calculate_next_occurrence)
        self.depleted = []
        super(FirstReactionMethod, self).__init__(process, state, **sampler_options)

    @property
    def transitions(self):
        return self.firings.keys()

    def update_state(self, dct):
        super(NextReactionMethod, self).update_state(dct)
        affected = self.dependency_graph.affected_transitions(dct)
        self.firings.update_items(affected)

    def add_transition(self, transition, **params):
        self.dependency_graph.add_reaction(transition)
        self.firings.add_item(transition, **params)

    def prune_transitions(self):
        if not self.depleted:
            return
        for trans in self.depleted:
            self.dependency_graph.remove_reaction(trans)
            self.firings.remove_item(trans)
        self.depleted = []

    def propose_potential_transition(self):
        if self.firings:
            return self.firings.topitem()
        else:
            return float('inf'), None, ()

    def perform_transition(self, time, transition, *args):
        super(NextReactionMethod, self).perform_transition(time, transition, *args)
        self.update_firing_times(time, transition, *args)

    def reject_transition(self, time, transition, *args):
        self.time = time
        transition.last_occurrence = time
        self.firings.update_one_instance(transition)

    def update_firing_times(self, time, transition, *args):
        """Update next occurrences of firing and all affected transitions"""
        # update affected firing times
        affected = self.dependency_graph.affected_transitions(transition.affected_species)
        self.firings.update_one_instance(transition)
        # XXX order in which reactions are updated should be preserved over independant runs
        self.firings.update_items(affected)
        # mark depleted reactions
        for trans in affected:
            if (any(occurrence[0] == float('inf')
                    for occurrence in self.firings[trans])
                    and (trans.rule or isinstance(trans, Event))):
                self.depleted.append(trans)

    def calculate_next_occurrence(self, transition, data):
        """Calculate next occurrence of given reaction."""
        return transition.next_occurrence(self.time, self.state, self.rng)


class AndersonMethod(NextReactionMethod):
    """Next reaction method modified for time-dependent processes

    A stochastic simulation algorithm, published in
    D. F. Anderson, J. Chem. Phys. 127, 214107 (2007)
    as Algorithm (3) 'modified next reaction method'.

    This sampler correctly treats non-autonomous transitions, i.e.
    transitions with time dependent stochastic rates.

    See help(StochasticSimulationAlgorithm) for usage information.
    """
    def add_transition(self, transition):
        """Add transition with own internal clock (T, P)"""
        from math import log
        super(AndersonMethod, self).add_transition(
            transition,
            T=0, P=-log(self.rng.random()), t=self.time
        )

    def perform_transition(self, time, transition, data):
        from math import log

        def int_a_dt(trans, delta_t):
            """Integrate propensity for given delta_t"""
            if isinstance(trans, Event):
                return 0
            else:
                return trans.propensity_integral(self.state, self.time, delta_t)

        data.T = data.P
        data.P -= log(self.rng.random())
        data.t = time

        affected = self.dependency_graph.affected_transitions(
            transition.affected_species
        )
        affected.discard(transition)

        for trans in affected:
            for t, data in self.firings[trans]:
                data.T += int_a_dt(trans, time-data.t)
                data.t = time

        super(AndersonMethod, self).perform_transition(time, transition)

    def calculate_next_occurrence(self, transition, data):
        """Determine next firing time of a transition in global time scale"""
        target = data.P-data.T
        if isinstance(transition, Event):
            return transition.next_occurrence(self.time)
        else:
            return self.time + transition.propensity_meets_target(
                self.state, self.time, target)

    def update_state(self, dct):
        # XXX AndersonMethod.update_state does not validate against DSMTS.
        warnings.warn("""AndersonMethod.update_state behavior is currently invalid.
        
        If this functionality is vital, use AndersonFRM instead,
        which provides a correct but significantly slower implementation.
        """, Warning)
        super(AndersonMethod, self).update_state(dct)


class AndersonFRM(FirstReactionMethod):
    """First reaction method modified for time-dependent processes

    This sampler is an inefficient implementation of Anderson's method,
    and will be removed in future versions.

    It should only be used when simulating processes with time-dependent
    transitions in situations where
    StochasticSimulationAlgorithm.update_state is used, because the
    more efficient implementation in AndersonMethod currently fails
    validation for this scenario.
    """
    def __init__(self, process, state, **sampler_options):
        warnings.warn("AndersonFRM will be removed in future versions.", DeprecationWarning)
        self.T = []
        self.P = []
        super(AndersonFRM, self).__init__(process, state, **sampler_options)

    def add_transition(self, transition):
        from math import log

        super(AndersonFRM, self).add_transition(transition)
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
        super(AndersonFRM, self).perform_transition(time, transition)


class CaoMethod(DirectMethod):
    """An approximate, explicit tau leap algorithm

    The algorithm has been published in

    Y. Cao, D. T. Gillespie, L. R. Petzold
    Efficient step size selection for the tau-leaping simulation method
    J. Chem. Phys. 124, 044109 (2006)

    Precision of the approximation can be controlled using the
    constructor keyword argument `epsilon`.

    See help(StochasticSimulationAlgorithm) for usage information

    Note: unlike other StochasticSimulationAlgorithm's iteration yields
    a tdictionary of ransition counts. 
    """
    n_crit = 10
    tauleap_threshold = 10
    micro_steps = 100

    def __init__(self, process, state, epsilon=0.03, tstart=0., tmax=float('inf'), steps=None, seed=None):
        from numpy.random import RandomState
        super(CaoMethod, self).__init__(process, state, tstart=tstart, tmax=tmax, steps=steps, seed=seed)
        self.num_reactions = 0
        self.epsilon = epsilon
        self.rng2 = RandomState(seed)
        self.abandon_tauleap = -1

    def __iter__(self):
        # XXX should become regular iterator
        from math import log

        while not self.has_reached_end():
            if self.time == self.tmax:
                break
            # step 1: partition reactions -- Eq. (10)
            Jcrit, Jncr = self.identify_critical_reactions()

            Irs = set(s for trans in self.transitions for s in trans.reactants)
            Incr = set(s for trans in Jncr for s in trans.reactants)

            if Incr:
                # step 2: determine noncritical tau -- Eqs. (32) and (33)
                mu = {s: sum(trans.stoichiometry.get(s, 0)*a
                             for trans, a in Jncr.items()) for s in Incr}
                var = {s: sum(trans.stoichiometry.get(s, 0)**2*a
                              for trans, a in Jncr.items()) for s in Incr}
                eps = {s: max(self.epsilon*self.state[s]*self.gi(s), 1.) for s in Incr}

                tau_ncr1 = min((eps[s]/abs(mu[s])) if mu[s] else float('inf') for s in Incr)
                tau_ncr2 = min((eps[s]**2/abs(var[s])) if var[s] else float('inf') for s in Incr)
                tau_ncr = min((tau_ncr1, tau_ncr2, self.tmax-self.time))
            else:
                tau_ncr = float('inf')

            a0 = sum(mult*prop for _, prop, mult in self.propensities.items())

            if not a0:
                break

            while True:
                # step 3: abandon tau leaping if not enough expected gain
                if tau_ncr <= self.tauleap_threshold / a0:
                    if self.abandon_tauleap == -1:
                        self.abandon_tauleap = self.step
                    it = DirectMethod.__iter__(self)
                    for _ in range(self.micro_steps):
                        try:
                            trans = next(it)
                        except StopIteration:
                            return
                        self.num_reactions += 1
                        yield trans
                    break
                elif self.abandon_tauleap != -1:
                    logging.debug("Abandoned tau-leaping for %d steps" % (self.step-self.abandon_tauleap))
                    self.abandon_tauleap = -1

                # step 4: determine critical tau
                ac = sum(propensity for trans, propensity in Jcrit.items())
                tau_crit = -log(self.rng.random())/ac if ac else float('inf')

                # step 5: determine actual tau
                tau = min((tau_ncr, tau_crit, self.tmax-self.time))

                # step 5a
                firings = {trans: self.rng2.poisson(propensity*tau)
                           for trans, propensity in Jncr.items()}
                firings = {trans: n for trans, n in firings.items() if n}

                # step 5b
                if tau == tau_crit:
                    # fire exactly one ciritical reaction
                    transition = None
                    pick = self.rng.random()*ac
                    for transition, propensity in Jcrit.items():
                        pick -= propensity
                        if pick < 0.:
                            break
                    firings[transition] = 1

                new_reactions = sum(firings.values())

                # avoid overshooting self.steps
                if self.steps and self.step+new_reactions > self.steps:
                    tau_ncr /= 2
                    continue

                all_reactants = sum((n*trans.true_reactants
                                     for trans, n in firings.items()), multiset())
                all_products = sum((n*trans.true_products
                                    for trans, n in firings.items()), multiset())
                net_reactants = all_reactants - all_products
                net_products = all_products - all_reactants

                # step 6a: avoid negative copy numbers
                if any(self.state[s]<n
                       for s, n in net_reactants.items()):
                    tau_ncr /= 2
                    continue

                else:
                    # step 6b: perform transitions
                    self.time += tau
                    self.step += 1
                    self.num_reactions += new_reactions
                    self.state -= net_reactants
                    for rule in self.process.rules:
                        for trans in rule.infer_transitions(net_products, self.state):
                            trans.rule = rule
                            self.add_transition(trans)
                    self.state += net_products

                    # update propensities
                    affected_species = set().union(*(trans.affected_species
                                                     for trans in firings))
                    affected = self.dependency_graph.affected_transitions(affected_species)
                    self.update_propensities(affected)
                    self.prune_transitions()

                    yield None # XXX should yield AggregateTransition
                    break

        if self.step != self.steps and self.tmax < float('inf'):
            self.time = self.tmax

    def identify_critical_reactions(self):
        def critical(trans, propensity):
            if not propensity:
                return False
            elif not trans.true_reactants:
                return False
            elif L(trans) < self.n_crit:
                return True
            elif any(self.state[s] < rule.order for rule in self.process.rules for s in trans.true_products):
                return True
            else:
                return False

        def L(trans):
            return min(-float(self.state[s])/trans.stoichiometry.get(s, 0)
                       for s in trans.true_reactants
                       if trans.stoichiometry.get(s, 0) < 0)

        crit = {}
        noncrit = {}
        for trans, a, mult in self.propensities.items():
            partition = crit if critical(trans, a) else noncrit
            partition[trans] = mult*a
        return crit, noncrit

    def gi(self, species):
        g = 0
        for trans in self.transitions:
            if species not in trans.reactants:
                continue
            order = sum(trans.reactants.values())
            if order == 1:
                g = max(g, 1)
            elif order == 2 and trans.reactants[species] == 1:
                g = max(g, 2)
            elif order == 2 and trans.reactants[species] == 2:
                g = max(g, 2+(self.state[species]-1)**-1) if self.state[species] >= 2 else 3
            elif order == 3 and trans.reactants[species] == 1:
                g = max(g, 3)
            elif order == 3 and trans.reactants[species] == 2:
                g = max(g, 1.5*(2+(self.state[species]-1)**-1)) if self.state[species] >= 2 else 4.5
            elif order == 3 and trans.reactants[species] == 3:
                g = max(g, 3+(self.state[species]-1)**-1+(self.state[species]-2)**-1) if self.state[species] >= 3 else 5.5
            else:
                raise RuntimeError("Tau-leaping not implemented for reactions of order %d" % order)
        return g


if sys.argv[0].endswith('stocal/examples/validation.py'):
    # Provide some specialized CaoMethod samplers with fixed epsilon values

    class CaoMethod_003(CaoMethod):
        def __init__(self, process, state, tstart=0., tmax=float('inf'), steps=None, seed=None):
            super(CaoMethod_003, self).__init__(process, state, 0.03, tstart=tstart, tmax=tmax, steps=steps, seed=seed)
    
    class CaoMethod_001(CaoMethod):
        def __init__(self, process, state, tstart=0., tmax=float('inf'), steps=None, seed=None):
            super(CaoMethod_001, self).__init__(process, state, 0.01, tstart=tstart, tmax=tmax, steps=steps, seed=seed)
    
    class CaoMethod_0003(CaoMethod):
        def __init__(self, process, state, tstart=0., tmax=float('inf'), steps=None, seed=None):
            super(CaoMethod_0003, self).__init__(process, state, 0.003, tstart=tstart, tmax=tmax, steps=steps, seed=seed)
    
    class CaoMethod_0001(CaoMethod):
        def __init__(self, process, state, tstart=0., tmax=float('inf'), steps=None, seed=None):
            super(CaoMethod_0001, self).__init__(process, state, 0.001, tstart=tstart, tmax=tmax, steps=steps, seed=seed)
    
    class CaoMethod_00003(CaoMethod):
        def __init__(self, process, state, tstart=0., tmax=float('inf'), steps=None, seed=None):
            super(CaoMethod_00003, self).__init__(process, state, 0.0003, tstart=tstart, tmax=tmax, steps=steps, seed=seed)
    
    class CaoMethod_00001(CaoMethod):
        def __init__(self, process, state, tstart=0., tmax=float('inf'), steps=None, seed=None):
            super(CaoMethod_00001, self).__init__(process, state, 0.0001, tstart=tstart, tmax=tmax, steps=steps, seed=seed)
