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

import abc
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
        warnings.warn("Method will return an iterator in future versions.", DeprecationWarning)
        return [
            key for key, (prop, mult) in self._dict.items()
            for i in range(mult)
        ]

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
        warnings.warn("Method will return an iterator in future versions.", DeprecationWarning)
        return [
            key for key, instances in self._queue.items()
            for i in instances
        ]

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

        self.state = multiset()
        self.update_state(state)

        for transition in self.process.transitions:
            self.add_transition(transition)

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

    @abc.abstractproperty
    def transitions(self):
        """Return list of all transitions

        Implementations need to return every instance of an added
        transition, i.e. two copies of a transition that has
        multiplicity two.
        """
        return []

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
        by TrajectorySampler.propose_potential_transition.
        """
        self.step += 1
        self.time = time # time of transition, in order to perform transition
        transition.last_occurrence = time
        self.state -= transition.true_reactants # Takes away the used up reactant from their respective total species count
        for rule in self.process.rules:
            for trans in rule.infer_transitions(transition.true_products, self.state):
                trans.rule = rule
                self.add_transition(trans)
        self.state += transition.true_products # Adds products to their respective total species

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
                self.prune_transitions()
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
        self.dependency_graph = DependencyGraph()
        self.propensities = MultiDict()
        self.depleted = []
        super(DirectMethod, self).__init__(process, state, t, tmax, steps, seed)

    def add_transition(self, transition): # not loop related
        self.dependency_graph.add_reaction(transition) # calls method(in this case the method/function is "dependency_graph.add_reaction()")  of inherited class(superclass) for any given object of this class(represented as self)
        propensity = self.calculate_propensity(transition) # calls method(in this case the method/function is "calculate_propensity()")  of inherited class(superclass) for any given object of this class(represented as self)
        self.propensities.add_item(transition, propensity) # calls method(in this case the method/function is "propensities.add_item()")  of inherited class(superclass) for any given object of this class(represented as self)

    def update_state(self, dct): # Deals with concentrations of chemical species in system.
        super(DirectMethod, self).update_state(dct)
        affected_transitions = self.dependency_graph.affected_transitions(dct)
        self.update_propensities(affected_transitions)

    def prune_transitions(self): # Deals with concentrations of chemical species in system. Takes away reaction with propensities = 0
        for trans in self.depleted:
            self.dependency_graph.remove_reaction(trans) # removes transition from dependency graph
            del self.propensities[trans] # Deleting both transition and it's propensities in the dictionary(MultiDict)
        self.depleted = []

    @property
    def transitions(self):
        return self.propensities.keys()

    def propose_potential_transition(self): # relates loop = pseudocode
        from math import log

        total_propensity = sum(mult*prop
                               for transition, prop, mult
                               in self.propensities.items()) # relates to step 6
        if not total_propensity:
            return float('inf'), None, tuple()

        delta_t = -log(self.rng.random())/total_propensity # step 2 of the linear time algorithm

        transition = None
        pick = self.rng.random()*total_propensity
        for transition, prop, mult in self.propensities.items():
            pick -= mult*prop # pick = pick - mult*prop
            if pick < 0.:
                break

        return self.time + delta_t, transition, tuple()

    def perform_transition(self, time, transition): # relates to step 4.
        super(DirectMethod, self).perform_transition(time, transition)
        affected = self.dependency_graph.affected_transitions(transition.affected_species)
        self.update_propensities(affected)

    def update_propensities(self, affected_transitions): # relates to step 4 part 1.
        """Update propensities of all given transitions.

        If the new propensity of a transition is 0 and the transition
        has been derived by a rule or is an Event, the transition gets
        added to self.depleted for later pruning."""
        for trans in affected_transitions:
            propensity = trans.propensity(self.state)
            if propensity == 0 and (trans.rule or isinstance(trans, Event)):
                self.depleted.append(trans)
            self.propensities.update_item(trans, propensity)

    def calculate_propensity(self, transition):
        """Return propensity of the given transition for current state."""
        return transition.propensity(self.state)


class FirstReactionMethod(TrajectorySampler):
    """Implementation of Gillespie's first reaction method.

    A stochastic simulation algorithm, published in
    D. T. Gillespie, J. Comp. Phys. 22, 403-434 (1976).

    FirstReactionMethod works with processes that feature deterministic
    transitions, i.e. Event's.

    See help(TrajectorySampler) for usage information.
    """
    def __init__(self, process, state, t=0., tmax=float('inf'), steps=None, seed=None):
        self._transitions = []
        self.firings = []
        super(FirstReactionMethod, self).__init__(process, state, t, tmax, steps, seed)

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
            return transition.reactants <= self.state # isn't applicable
        else:
            return super(FirstReactionMethod, self).is_applicable(time, transition, *args) # is applicable

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

    See help(TrajectorySampler) for usage information.
    """
    def __init__(self, process, state, t=0., tmax=float('inf'), steps=None, seed=None):
        self.dependency_graph = DependencyGraph()
        self.firings = PriorityQueue(self.calculate_next_occurrence)
        self.depleted = []
        super(FirstReactionMethod, self).__init__(process, state, t, tmax, steps, seed)

    def add_transition(self, transition, **params):
        self.dependency_graph.add_reaction(transition)
        self.firings.add_item(transition, **params)

    def update_state(self, dct):
        super(NextReactionMethod, self).update_state(dct)
        affected = self.dependency_graph.affected_transitions(dct)
        self.firings.update_items(affected)

    def prune_transitions(self):
        for trans in self.depleted:
            self.dependency_graph.remove_reaction(trans)
            self.firings.remove_item(trans)
        self.depleted = []

    @property
    def transitions(self):
        return self.firings.keys()

    def propose_potential_transition(self):
        if self.firings:
            return self.firings.topitem()
        else:
            return float('inf'), None, ()

    def perform_transition(self, time, transition, *args):
        super(NextReactionMethod, self).perform_transition(time, transition, *args)
        self.update_firing_times(time, transition, *args)

    def update_firing_times(self, time, transition, *args):
        """Update next occurrences of firing and all affected transitions"""
        # update affected firing times
        affected = self.dependency_graph.affected_transitions(transition.affected_species)
        self.firings.update_one_instance(transition)
        self.firings.update_items(affected)
        # mark depleted reactions
        for trans in affected:
            if (any(occurrence[0] == float('inf')
                    for occurrence in self.firings[trans])
                    and (trans.rule or isinstance(trans, Event))):
                self.depleted.append(trans)

    def reject_transition(self, time, transition, *args):
        self.time = time
        transition.last_occurrence = time
        self.firings.update_one_instance(transition)

    def calculate_next_occurrence(self, transition, data):
        """Calculate next occurrence of given reaction."""
        return transition.next_occurrence(self.time, self.state, self.rng)


class AndersonNRM(NextReactionMethod):
    """Next reaction method modified for time-dependent processes

    A stochastic simulation algorithm, published in
    D. F. Anderson, J. Chem. Phys. 127, 214107 (2007)
    as Algorithm (3) 'modified next reaction method'.

    This sampler correctly treats non-autonomous transitions, i.e.
    transitions with time dependent stochastic rates.

    See help(TrajectorySampler) for usage information.
    """
    def add_transition(self, transition):
        """Add transition with own internal clock (T, P)"""
        from math import log
        super(AndersonNRM, self).add_transition(
            transition,
            T=0, P=-log(self.rng.random())
        )

    def calculate_next_occurrence(self, transition, data):
        """Determine next firing time of a transition in global time scale"""
        target = data.P-data.T
        if isinstance(transition, Event):
            return transition.next_occurrence(self.time)
        else:
            return self.time + transition.propensity_meets_target(
                self.state, self.time, target)

    def perform_transition(self, time, transition, data):
        from math import log

        def int_a_dt(trans, delta_t):
            """Integrate propensity for given delta_t"""
            if isinstance(trans, Event):
                return 0
            else:
                return trans.propensity_integral(self.state, self.time, delta_t)

        data.P -= log(self.rng.random())
        affected = self.dependency_graph.affected_transitions(
            transition.affected_species
        )
        for trans in affected:
            for time, data in self.firings[trans]:
                data.T += int_a_dt(trans, time-self.time)

        super(AndersonNRM, self).perform_transition(time, transition)

class Group(object): # holds information for each group of transitions

    def __init__(self, gmin, gmax):
        self.gmin = gmin
        self.gmax = gmax
        self.num_transitions = 0
        self.transitions = []
        self.propensities = []

    def group_total_propensities(self): # sum of all propensities for the transition in the group
        return sum(self.propensities)

    def transition_propensity(self, r): # returns the propensity for a given transition
        l = self.transitions.index(r)
        propensity = self.propensities[l]
        return propensity

class CompositionRejection(TrajectorySampler):
    """Implementation of Composition and Rejection method."""
    PMIN = 0.03125 # 1/32 = 0.03125
    PMAX = 1024 # (1/32) x 2^15 = 1024

    def __init__(self, process, state, t=0., tmax=float('inf'), steps=None, seed=None):
        if any(not isinstance(r, Reaction) for r in process.transitions):
            raise ValueError("CompositionRejection only works with Reactions.")
        if any(not issubclass(r.Transition, Reaction) for r in process.rules):
            raise ValueError("CompositionRejection only works with Reactions.")
        self.dependency_graph = DependencyGraph()
        self.depleted = []
        self.groups = []
        self.gzero = Group(0, 0)
        self.gmin = self.PMIN
        self.gmax = 2*self.PMIN
        self.group_index = {} # empty dictionary, holding group and index of pruned(self.depleted) reactions
        super(CompositionRejection, self).__init__(process, state, t, tmax, steps, seed)

        while self.gmax < self.PMAX: # this while loop groups the reactions into groups
            new_group = Group(self.gmin, self.gmax) # a object of the class Group is created and called "new_group"
            self.groups.append(new_group) # self.groups is a list, which has object called "new_group" added as it's item
            self.gmin *= 2 # gmin = gmin * 2, for every iteration of loop gmin value doubles
            self.gmax *= 2 # gmax = gmax * 2, for every iteration of loop gmax value doubles


    def add_transition(self, transition): # adds transitions
        prop = self.calculate_propensity(transition) # assigns propensity to prop variable
        self.add_transition_with_prop(transition, prop)
        self.dependency_graph.add_reaction(transition) # adds transition to dependency graph


    def add_transition_with_prop(self, transition, propensity): # Sorts transitions into either Gzero or Self.groups
        prop = propensity
        assert self.PMIN <= prop <= self.PMAX or prop == 0, "Abnormal prop"

        if prop == 0:
            g = self.gzero

            g.transitions.append(transition)
            g.propensities.append(prop)
            g.num_transitions += 1
        else:
            for grp in self.groups:
                if grp.gmax >= prop > grp.gmin:
                    gr = grp

                    gr.transitions.append(transition)
                    gr.propensities.append(prop)
                    gr.num_transitions += 1
                    break # terminates the current loop and resumes execution at the next statement.


    def update_state(self, dct): # Deals with concentrations of chemical species in system.
        super(CompositionRejection, self).update_state(dct)
        affected_transitions = self.dependency_graph.affected_transitions(dct)
        self.update_propensities(affected_transitions)

    def prune_transitions(self): # Takes away reaction with propensities of 0.
        for trans in self.depleted:
            self.dependency_graph.remove_reaction(trans)
            del self.propensities[trans]
            g = self.group_index[trans]
            Group = g[0] # the last group the transition was in
            index = g[1] # the list index of transition and propensity in Group
            del Group.propensities[index] # removes trans's propensity from Group propensity list
            del Group.transitions[index] # removes trans from Group transition list

        self.depleted = [] # empties given list
        self.group_index = {} # empties given list

    @property
    def transitions(self): # outputs list of all current transitions
        a = [x.transitions for x in self.groups] # list of all transitions in self.groups
        b = self.gzero.transitions # self.gzero's transition list assigned to variable b.
        q = [x for x in a if x not in b] # list of transitions in self.groups that aren't in self.gzero
        all_transitions = b + q # Transitions in both self.groups and self.gzero
        return all_transitions

    def propose_potential_transition(self):  # proposes a transition, relates to Composition aspect of SSA-CR
        from math import log

        total_propensity = sum(g.group_total_propensities() for g in self.groups)

        if not total_propensity:
            return float('inf'), None, (None, )

        import random
        pick = random.uniform(0, 1) * total_propensity # each reaction has an equal chance(likelihood) of occurring, "r1"
        for group in self.groups:
            pick -= group.num_transitions
            if pick <= 0:
                group_selected = group # the chosen group
                break

        delta_t = -log(self.rng.random()) / total_propensity # computes time increment, "r2"

        r3 = random.randint(1, group_selected.num_transitions) # x-axis

        transition = group_selected.propensities[r3-1]  # picks a reaction/transition



        return self.time + delta_t, transition, (group_selected, )

    def perform_transition(self, time, transition, group_selected): # relates to step 4 of SSA-CR algorithm.
        if transition == None:
            pass
        else:
            super(CompositionRejection, self).perform_transition(time, transition, group_selected)
            affected = self.dependency_graph.affected_transitions(transition.affected_species)
            self.update_propensities(affected)



    def is_applicable(self, time, transition, group_selected): # relates to the rejection aspect of SSA-CR algorithm

        if transition == None:
            pass
        else:
            m = group_selected.gmax
            propensity = group_selected.transition_propensity(transition)
            import random
            r4 = random.uniform(0, m) # y-axis, computes r4
            if propensity >= r4: # is_applicable returns True, if r4 =< propensity
                return False  # isn't applicable
            else:
                return True # is applicable


    def update_propensities(self, affected_transitions): # relates to step 4 part 1. Deals with regrouping.
        """Update propensities of all given transitions.

        If the new propensity of a transition is 0 and the transition
        has been derived by a rule or is an Event, the transition gets
        added to self.depleted for later pruning."""

        for g in self.groups: # Looks at all affected reactions and updates their propensity
            for trans in affected_transitions:
                if trans in g.transitions:
                    prop = trans.propensity(self.state) # assigns update propensity to prop variable
                    t = g.propensities.index(g.transition_propensity(trans)) # Find's the index of the transition's propensity
                    g.propensities[t] = prop  # Updates the propensity of the transition in it's list
                    if g.gmin < prop <= g.max: # build on <--
                        pass  # update trans propensity value in the group
                    else: # # if transition doesn't belong to it's group due to prop change, it is regrouped

                        g.transitions.remove(trans) # removes transitions from it's previous group
                        g.propensities.remove(prop) # removes transition's propensity from it's previous group
                        g.num_transitions -= 1 # reduces number of group transitions by one
                        self.add_transition_with_prop(trans, prop) # performs the regrouping


        for trans in affected_transitions:
            if trans in self.gzero.transitions:
                prop = trans.propensity(self.state) # assigns update propensity to prop variable
                t = self.gzero.propensities.index(self.gzero.transition_propensity(trans))  # Find's the index of the transition's propensity
                self.gzero.propensities[t] = prop  # Updates the propensity of the transition in it's list
                if self.gzero.gmin < prop <= self.gzero.gmax: # build on <--
                    pass  # update trans' propensity value in the group
                else: # # if transition doesn't belong to it's group due to prop change, it is regrouped

                    self.gzero.transitions.remove(trans) # removes transitions from it's previous group
                    self.gzero.propensities.remove(prop) # removes transition's propensity from it's previous group
                    self.gzero.num_transitions -= 1 # reduces number of group transitions by one
                    self.add_transition_with_prop(trans, prop) # performs the regrouping

        for trans in affected_transitions: # prunes affected reaction if prop = 0
            propensity = trans.propensity(self.state)
            if propensity == 0:
                for g in self.groups:
                    if trans in g.transitions:
                        t = g.transitions.index(trans)
                        index = g.propensities[t]
                        group = g
                        self.group_index[trans] = [group, index]
                        self.depleted.append(trans)

    def calculate_propensity(self, transition):
        """Return propensity of the given transition for current state."""
        return transition.propensity(self.state)






