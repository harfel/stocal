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
"""Transitions and rules for the stocal framework

This module provides classes that model transitions of chemicals,
such as reactions and events, as well as rules that can derive
novel transitions during the course of a stochastic simulation.
"""

import abc
import warnings

from .utils import with_metaclass
from .structures import multiset


class Transition(with_metaclass(abc.ABCMeta, object)):
    """Abstract base class for transformations of reactants into products.

    Transitions define reactants and products, but not the kinetic
    law of the transformation. Consult MassAction and Event for the most
    common rate law implementations.

    Transition instances provide the attributes self.reactants and
    self.products, which are mappings of chemical species to their
    stoichiometric factors. Instances also provide the attributes
    self.true_reactants and self.true_products, which exclude chemical
    species that appear both as reactants and products (i.e. as
    catalysts).

    Instances also have an attribute self.rule which defaults to None
    and which is set by the TrajectorySampler if a Transition has been
    inferred by a Rule.

    Modiying any of these attributes after initialization is an error
    and leads to undefined behavior.
    """

    rule = None

    def __init__(self, reactants, products):
        """Initialization

        reactants and products are either mappings that give the
        stoichiometric factor of each involved species, or sequences
        which are interpreted as unordered lists.
        """
        reactants = multiset(reactants)
        products = multiset(products)

        if not all(n > 0 for n in reactants.values()):
            raise ValueError("reactant stoichiometries must be positive.")
        if not all(n > 0 for n in products.values()):
            raise ValueError("product stoichiometries must be positive.")
        if not reactants and not products:
            raise ValueError(
                "%s must have either reactants or products."
                % type(self).__name__
            )

        self.reactants = reactants
        self.products = products

        self.true_reactants = reactants - products
        self.true_products = products - reactants

        self.last_occurrence = -1.
        self._hash = 0

    def __eq__(self, other):
        """Structural congruence

        Transitions are equal if their reactants and products are equal.
        """
        return (
            isinstance(other, Transition) and
            self.reactants == other.reactants and
            self.products == other.products
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        if not self._hash:
            self._hash = hash((
                tuple(sorted(self.reactants.items())),
                tuple(sorted(self.products.items())),
            ))
        return self._hash

    def __repr__(self):
        try:
            return '%s(%s, %s)' % (
                type(self).__name__, self.reactants, self.products
            )
        except AttributeError:
            return super(Transition, self).__repr__()

    def __str__(self):
        def dct2str(dct):
            """render multiset as sum of elements"""
            return ' + '.join(
                str(s) if n == 1 else '%s*%s' % (n, s) for s, n in dct.items()
            )
        try:
            return '%s --> %s' % (dct2str(self.reactants), dct2str(self.products))
        except AttributeError:
            return super(Transition, self).__str__()

    @abc.abstractmethod
    def next_occurrence(self, time, state):
        """Time of next occurrence after given time.

        This method has to be implemented by a subclass. See
        Event.next_occurrence and Reaction.next_occurrence for examples.
        """
        return float('inf')


class Reaction(Transition):
    """Abstract base class for stochastic transitions.

    Stochastic transitions are those that occur with a certain
    proponsity within a given time interval.
    Subclasses must implement the propensity method. For autonomous
    reactions, i.e. where the propensity does only depend on the state
    but not on time, the propensity method must have the signature

    def propensity(self, state):
        ...

    If the propensity of a reaction does depend on time, the propensity
    method must have the signature

    def propensity(self, state, time):
        ...

    In the latter case, the user can additionally override the methods
    Reaction.propensity_integral and Reaction.propensity_meets_target,
    for example with analytic solutions.
    """
    def __eq__(self, other):
        """Structural congruence

        Reactions are equal if their reactants, products, and propensity
        functions are equal.
        """
        return (
            super(Reaction, self).__eq__(other) and
            isinstance(other, Reaction) and
            type(self).propensity == type(other).propensity
        )

    def __hash__(self):
        if not self._hash:
            self._hash = hash((
                super(Reaction, self).__hash__(),
                type(self).propensity
            ))
        return self._hash

    @property
    def is_autonomous(self):
        """True if propensity does not depend on time.

        The value of this property is inferred from the signature
        of the propensity method. See help(Reaction) for details."""
        import inspect
        try:
            # python 3
            return len(inspect.signature(self.propensity).parameters) == 1
        except AttributeError:
            # python 2.7
            return len(inspect.getargspec(self.propensity).args) == 2

    @abc.abstractmethod
    def propensity(self, state):
        """Reaction propensity.

        This method has to be provided by a subclass and must return a
        non-negative float denoting the reaction propensity for the
        provided state. The function should not modify the provided state.
        """
        return 0.

    def next_occurrence(self, time, state):
        """Determine next reaction firing time.

        This is a helper function to use Reactions in next-firing-time
        based TrajectorySampler's. The method randomly draws a delay
        from a Poisson distribution with mean propensity  and returns
        the given current time plus the delay.
        """
        from random import random
        from math import log
        propensity = self.propensity(state)
        if not propensity:
            return float('inf')
        else:
            return time - log(random())/propensity

    def propensity_integral(self, state, time, delta_t):
        """Integrate propensity function from time to time+delta_t

        If the reaction is autonomous, this method returns
        self.propensity(state)*delta_t. If it is non-autonomous,
        the integral is calculated numerically. Override this method
        if the propensity function can be integrated analytically.
        (Requires scipy).
        """
        if self.is_autonomous:
            return self.propensity(state)*delta_t
        elif delta_t == float('inf'):
            return float('inf')
        else:
            from scipy.integrate import quad as integral
            propensity = lambda t: self.propensity(state, t)
            return integral(propensity, time, time+delta_t)[0]

    def propensity_meets_target(self, state, time, target):
        """Time window at which propensity integral meets a given target

        If the reaction is autonomous, this method returns
        target/self.propensity(state). If it is non-autonomous, the
        time window dt is numerically evaluated as the interval dt
        at which

        integral_t^{t+dt} a(X(t), s) ds = target.

        Override this method if an analytic solution to the integral
        can be obtained.
        (Requires scipy).
        """
        if self.is_autonomous:
            propensity = self.propensity(state)
            return target/propensity if propensity else float('inf')
        elif self.reactants not in state:
            return float('inf')
        else:
            from scipy.optimize import minimize_scalar as minimize
            fun = lambda dt: (self.propensity_integral(state, time, dt)-target)**2
            return minimize(fun).x


class MassAction(Reaction):
    """Reactions with mass action kinetics.

    The propensity of a mass action reaction is defined as the
    stochastic rate constant c times the number of possible reaction
    partners, where the latter is the product of the binomial
    coefficients to choose n reaction partners out of m molecules for
    each reacting species.
    """
    def __init__(self, reactants, products, c):
        if c < 0:
            raise ValueError("stochastic rate constants must be non-negative.")
        super(MassAction, self).__init__(reactants, products)
        self.constant = c

    def __repr__(self):
        try:
            return '%s(%s, %s, %g)' % (
                type(self).__name__, self.reactants, self.products, self.c
            )
        except AttributeError:
            return super(MassAction, self).__repr__()

    def propensity(self, state):
        """Reaction propensity for the given state.

        Calling propensity does not modify the underlying reaction.
        """
        from functools import reduce # for python3 compatibility

        if not isinstance(state, multiset):
            warnings.warn("state must be a multiset.", DeprecationWarning)

        def choose(n, k):
            """binomial coefficient"""
            return reduce(lambda x, i: x*(n+1-i)/i, range(1, k+1), 1)
        return reduce(
            lambda a, b: a*b,
            (choose(state.get(s, 0), n) for s, n in self.reactants.items()),
            self.constant
        )

    def __eq__(self, other):
        """Structural congruence

        MassAction reactions are equal if their reactants, products,
        and stochastic rate constants are equal.
        """
        return (
            super(MassAction, self).__eq__(other) and
            isinstance(other, MassAction) and
            self.constant == other.constant
        )

    def __hash__(self):
        if not self._hash:
            self._hash = hash((
                super(MassAction, self).__hash__(),
                self.constant
            ))
        return self._hash


class Event(Transition):
    """Deterministic transitions.

    Events are Transition's that occur either once at a specified time,
    or periodically with a given frequency starting at a specified time.
    """
    is_autonomous = True

    def __init__(self, reactants, products, time, frequency=0):
        if time < 0:
            raise ValueError("time must be greater than 0.")
        if frequency < 0:
            raise ValueError("dt must be greater than (or equal to) 0.")
        super(Event, self).__init__(reactants, products)
        self.time = time
        self.frequency = frequency

    def __eq__(self, other):
        """Structural congruence

        Events are equal if their reactants, products, time, and
        frequency are equal.
        """
        return (
            super(Event, self).__eq__(other) and
            isinstance(other, Event) and
            self.time == other.time and
            self.frequency == other.frequency
        )

    def __hash__(self):
        if not self._hash:
            self._hash = hash((
                super(Event, self).__hash__(),
                self.time, self.frequency
            ))
        return self._hash

    def next_occurrence(self, time, state=None):
        """Next occurrence of the Event at or after time.

        If the event does not re-occur, returns float('inf').
        Calling next_occurrence leaves the event unmodified.
        """
        if time < self.time:
            return self.time
        elif self.frequency:
            future = time + (self.time-time)%self.frequency
            return (future
                    if self.last_occurrence != time
                    else future+self.frequency)
        elif time == self.time and self.last_occurrence != time:
            return time
        else:
            return float('inf')


class Rule(with_metaclass(abc.ABCMeta, object)):
    """Abstract base class for rules

    Subclasses must provide a class attribute called Transition,
    which denotes the class of Transitions they generate.
    They also must provide a method infer_transitions which performs
    the actual transition inference.
    """

    @abc.abstractproperty
    def Transition(self):
        """The type of Transition that the rule generates"""
        return Transition

    @abc.abstractmethod
    def infer_transitions(self, new_species, state):
        """infer new transitions among new_species and state

        This method is called by the stochastic simulation algorithm
        "in the middle of a transition", i.e. after removing reactants
        from the state but before adding new_species as products.
        Implementations must return an iterable of Transition objects.
        """
        # XXX Should this take last_transition?
        raise StopIteration


class ReactionRule(Rule):
    """Abstract base class that facilitates inference of Reactions

    This class provides a standard implementation of infer_transitions
    that generates all possible reactant combinations from species in
    state and new_species, that only became possible because of species
    in new_species, and could not have been formed by reactants in state
    alone.
    If ReactionRule.signature is given, it must evaluate to a sequence
    of type objects. Combinations are then only formed among reactants
    that are instances of the given type.
    For each combination, the inference algorithm calls
    ReactionRule.novel_reactions. This method, to be implemented by a
    subclass, should return an iterable over every reaction that takes
    the novel species as reactants.
    """

    @abc.abstractmethod
    def novel_reactions(self, *reactants):
        """Infer reactions for the given unordered list of reactants.

        To be implemented by a subclass.
        """
        raise StopIteration

    @property
    def Transition(self):
        """Return transition type from novel_reactions return annotation

        In python3, ReactionRule.Transition is optional and can alternatively
        be provided as novel_reactions return type annotation:

        from typing import Iterator
        class MyRule(stocal.ReactionRule):
            def novel_reactions(self, *reactants) -> Iterator[TransitionClass]:

        In python2, the property raises an AttributeError.
        """
        import inspect
        try:
            # python 3
            signature = inspect.signature(self.novel_reactions)
            ret_ann = signature.return_annotation
            cls = ret_ann.__args__[0]
            if not issubclass(cls, Transition):
                raise TypeError("%s is not a subclass of stocal.Transition"
                                % cls.__name__)
            else:
                return cls
        except AttributeError:
            raise TypeError("%s.Transition not defined and not inferable"
                            +" from novel_reactions signature"
                            % type(self).__name__)

    @property
    def order(self):
        """Reaction order of infered reactions.

        The order of a reaction is the number of reactant molecules.
        To be defined by a subclass."""
        import inspect
        try:
            # python 3
            parameters = inspect.signature(self.novel_reactions).parameters
            return sum(1 for par in parameters.values()
                       if par.kind == inspect.Parameter.POSITIONAL_OR_KEYWORD)
        except AttributeError:
            # python 2.7
            return len(inspect.getargspec(self.novel_reactions).args) - 1

    @property
    def signature(self):
        """Type signature of ReactionRule.novel_reactions

        In python2, this defaults to self.order*[object]. Override the
        attribute in a subclass to constrain reactant types of
        novel_reactions.

        In python3, the signature is inferred from type annotations
        of the novel_reactions parameters (defaulting to object for
        every non-annotated parameter).
        """
        import inspect
        try:
            # python 3
            signature = inspect.signature(self.novel_reactions)
            return [p.annotation if p.annotation != inspect.Parameter.empty else object
                    for p in signature.parameters.values()]
        except AttributeError:
            return self.order*[object]

    def infer_transitions(self, new_species, state):
        """Standard inference algorithm for Reactions.

        see help(type(self)) for an explanation of the algorithm.
        """
        if not isinstance(new_species, multiset):
            warnings.warn("last_products must be a multiset.", DeprecationWarning)
        if not isinstance(state, multiset):
            warnings.warn("state must be a multiset.", DeprecationWarning)

        def combinations(reactants, signature, annotated_species, novel):
            """Yield all novel combinations comaptible with signature

            See class doc for details."""
            if not signature:
                if novel:
                    yield reactants
                return

            skipped = []
            while annotated_species:
                species, start, end = annotated_species.pop(0)
                if isinstance(species, signature[0]):
                    break
                else:
                    skipped.append((species, start, end))
            else:
                if not annotated_species:
                    return

            for combination in combinations(reactants,
                                            signature,
                                            skipped+annotated_species,
                                            novel):
                yield combination
            if end > 1:
                annotated_species.insert(0, (species, start-1, end-1))
            for combination in combinations(reactants+[species],
                                            signature[1:],
                                            skipped+annotated_species,
                                            novel or start == 1):
                yield combination

        # could be simplified if specification would enforce multiset state:
        # next_state = state + last_products
        # novel_species = sorted((
        #    (species, state[species]+1,
        #     min(next_state[species],
        #         len([typ for typ in self.signature if isinstance(species, typ)])))
        #    for species in set(new_species).union(state)
        #), key=lambda item: item[1]-item[2])

        novel_species = sorted((
            (species, state.get(species, 0)+1,
             min(new_species.get(species, 0)+state.get(species, 0),
                 len([typ for typ in self.signature if isinstance(species, typ)])))
            for species in set(new_species).union(state)
        ), key=lambda item: item[1]-item[2])
        for reactants in combinations([], self.signature, novel_species, False):
            for trans in self.novel_reactions(*reactants):
                yield trans



class Process(object):
    """Stochastic process class

    A collection of all transitions and rules that define a
    stochastic process. When initializing a TrajectorySampler with
    a Process instance, transitions get copied over into the sampler.
    This makes it possible to use a single process instance with
    multiple samplers.
    """
    def __init__(self, transitions=None, rules=None):
        self.transitions = transitions or []
        self.rules = rules or []

    def trajectory(self, state, t=0., tstart=0., tmax=float('inf'), steps=None):
        """Create trajectory sampler for given state

        The method automatically chooses a suitable sampler for the
        given stochastic process, initialized with the given state
        and time.
        """
        if t:
            warnings.warn("pass start time as tstart", DeprecationWarning)
        tstart = tstart or t

        def transition_types():
            """Yield all generated transtion types of the process"""
            for trans in self.transitions:
                yield trans
            for rule in self.rules:
                yield rule.Transition

        # DirectMethod for process with normal reactions
        if all(issubclass(r, Reaction) and r.is_autonomous
               for r in transition_types()):
            from .algorithms import DirectMethod as Sampler
        # FirstReactionMethod if all reactions are autonomous
        elif all(r.is_autonomous for r in transition_types()):
            from .algorithms import FirstReactionMethod as Sampler
        # AndersonNRM if reactions are non-autonomous
        else:
            from .algorithms import AndersonNRM as Sampler

        return Sampler(self, state, tstart, tmax, steps)
