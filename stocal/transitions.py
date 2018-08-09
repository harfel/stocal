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

from ._utils import with_metaclass
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
    and which is set by the StochasticSimulationAlgorithm if a
    Transition has been inferred by a Rule.

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
        self.affected_species = self.true_reactants.union(self.true_products).domain

        self.stoichiometry = {
            species: self.true_products[species]-self.true_reactants[species]
            for species in self.affected_species
        }

        self.last_occurrence = -float('inf')
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
            return '<%s %s>' % (
                type(self).__name__, self
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

    def next_occurrence(self, time, state, rng=None):
        """Determine next reaction firing time.

        This is a helper function to use Reactions in next-firing-time
        based StochasticSimulationAlgorithm's. The method randomly draws
        a delay from a Poisson distribution with mean propensity  and
        returns the given current time plus the delay.
        """
        from math import log
        if not rng:
            import random as rng
        propensity = self.propensity(state)
        if not propensity:
            return float('inf')
        else:
            return time - log(rng.random())/propensity

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
            return '<%s %s, %g>' % (
                type(self).__name__, self, self.constant
            )
        except AttributeError:
            return super(MassAction, self).__repr__()

    def propensity(self, state):
        """Reaction propensity for the given state.

        Calling propensity does not modify the underlying reaction.
        """
        from functools import reduce # for python3 compatibility

        state = state if isinstance(state, multiset) else multiset(state)

        def choose(n, k):
            """binomial coefficient"""
            return reduce(lambda x, i: x*(n+1-i)/i, range(1, k+1), 1)
        return reduce(
            lambda a, b: a*b,
            (choose(state[s], n) for s, n in self.reactants.items()),
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

    def __repr__(self):
        try:
            return '<%s %s, %g, frequency=%g>' % (
                type(self).__name__, self, self.time, self.frequency
            )
        except AttributeError:
            return super(Event, self).__repr__()

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

    def next_occurrence(self, time, state=None, rng=None):
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
        raise StopIteration


class _TransitionRuleMetaclass(abc.ABCMeta):
    def __call__(self):
        cls = super(_TransitionRuleMetaclass, self).__call__()
        cls.order = self.get_order(cls)
        if not hasattr(cls, 'signature'):
            cls.signature = self.get_signature(cls)
        if not hasattr(cls, 'Transition'):
            cls.Transition = self.get_Transition(cls)
        return cls

    def get_order(self, cls):
        """Reaction order of infered reactions.

        The order of a reaction is the number of reactant molecules.
        To be defined by a subclass."""
        import inspect
        try:
            # python 3
            parameters = inspect.signature(cls.novel_reactions).parameters
            return sum(1 for par in parameters.values()
                       if par.kind == inspect.Parameter.POSITIONAL_OR_KEYWORD)
        except AttributeError:
            # python 2.7
            return len(inspect.getargspec(cls.novel_reactions).args) - 1

    def get_signature(self, cls):
        """Type signature of TransitionRule.novel_reactions

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
            signature = inspect.signature(cls.novel_reactions)
            return [p.annotation if p.annotation != inspect.Parameter.empty else object
                    for p in signature.parameters.values()]
        except AttributeError:
            return cls.order*[object]

    def get_Transition(self, cls):
        """Return transition type from novel_reactions return annotation

        In python3, TransitionRule.Transition is optional and can
        alternatively be provided as novel_reactions return type annotation:

        from typing import Iterator
        class MyRule(stocal.TransitionRule):
            def novel_reactions(self, *reactants) -> Iterator[TransitionClass]:

        In python2, the property raises an AttributeError.
        """
        import inspect
        try:
            # python 3
            signature = inspect.signature(cls.novel_reactions)
            ret_ann = signature.return_annotation
            cls = ret_ann.__args__[0]
            if not issubclass(cls, Transition):
                raise TypeError("%s is not a subclass of stocal.Transition"
                                % cls.__name__)
            else:
                return cls
        except AttributeError:
            raise TypeError(("%s.Transition not defined and not inferable"
                            +" from novel_reactions signature")
                            % cls.__name__)


class TransitionRule(Rule, with_metaclass(_TransitionRuleMetaclass, Rule)):
    """Abstract base class that facilitates inference of Reactions

    This class provides a standard implementation of infer_transitions
    that generates all possible reactant combinations from species in
    state and new_species, that only became possible because of species
    in new_species, and could not have been formed by reactants in state
    alone. For each combination, the inference algorithm calls
    TransitionRule.novel_reactions. This method, to be implemented by a
    subclass, should return an iterable over every reaction that takes
    the novel species as reactants.

    If TransitionRule.signature is given, it must evaluate to a sequence
    of type objects. Combinations are then only formed among reactants
    that are instances of the given type. In python3, the signature can
    automatically be inferred from type annotations of the
    novel_reactions parameters (defaulting to object for every
    non-annotated parameter).
    
    In python3, TransitionRule.Transition can alternatively be provided
    as novel_reactions return type annotation:

        from typing import Iterator
        class MyRule(stocal.TransitionRule):
            def novel_reactions(self, *reactants) -> Iterator[TransitionClass]:
    """

    @abc.abstractmethod
    def novel_reactions(self, *reactants):
        """Infer reactions for the given unordered list of reactants.

        To be implemented by a subclass.
        """
        raise StopIteration

    def infer_transitions(self, new_species, state):
        """Standard inference algorithm for Reactions.

        see help(type(self)) for an explanation of the algorithm.
        """
        new_species = new_species if isinstance(new_species, multiset) else multiset(new_species)
        state = state if isinstance(state, multiset) else multiset(state)

        novel_species = sorted((
            (species, state[species]+1,
             min(new_species[species]+state[species],
                 len([typ for typ in self.signature if isinstance(species, typ)])))
            for species in set(new_species).union(state)
            if any(isinstance(species, typ) for typ in self.signature)
        ), key=lambda item: item[2]-item[1])
        for reactants in self._combinations([], self.signature, novel_species, False):
            for trans in self.novel_reactions(*reactants):
                yield trans

    def _combinations(self, reactants, signature, annotated_species, novel):
        """Yield all novel combinations comaptible with signature

        See class doc for details."""
        if not signature:
            if novel:
                yield reactants
            return
        elif not novel and all(end < start
                               for _, start, end in annotated_species):
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

        for combination in self._combinations(reactants,
                                              signature,
                                              skipped+annotated_species,
                                              novel):
            yield combination
        if end > 1:
            annotated_species.insert(0, (species, start-1, end-1))
        for combination in self._combinations(reactants+[species],
                                              signature[1:],
                                              skipped+annotated_species,
                                              novel or start == 1):
            yield combination


class Process(object):
    """Stochastic process class

    A collection of all transitions and rules that define a
    stochastic process. When initializing a StochasticSimulationAlgorithm
    with a Process instance, transitions get copied over into the sampler.
    This makes it possible to use a single process instance with
    multiple samplers.
    """
    def __init__(self, transitions=None, rules=None):
        self.transitions = transitions or []
        self.rules = rules or []

    def __eq__(self, other):
        return self.transitions == other.transitions and self.rules == other.rules

    def __ne__(self, other):
        return not self == other

    def trajectory(self, state, t=0., tstart=0., tmax=float('inf'), steps=None, seed=None):
        """Create trajectory sampler for given state

        Depreated: please use Process.sample instead.

        The method does the same as Process.trajectory, but returns
        a sampler that only yields individual Transition objects
        in each iteration:
 
        >>> process = Process()
        >>> state = {}
        >>> trajectory = process.trajectory(state)
        >>> for transition in trajectory():
        ...     print(trajectory.time, trajectory.state, transition)

        C.f. stocal.algorithms.StochasticSimulationAlgorithm for details
        on the sampler class.
        """
        warnings.warn("Use Process.sample instead", DeprecationWarning)
        return self._trajectory(state, tstart=tstart or t, tmax=tmax, steps=steps, seed=seed)

    def _trajectory(self, state, t=0., tstart=0., tmax=float('inf'), steps=None, seed=None):
        def transition_types():
            """Yield all generated transtion types of the process"""
            for trans in self.transitions:
                yield trans
            for rule in self.rules:
                yield rule.Transition

        # select suitable simulation algorithm
        if any(not r.is_autonomous for r in transition_types()):
            # AndersonMethod for processes with non-autonomous reactions
            from .algorithms import AndersonMethod as Sampler
        else:
            # NextReactionMethod for anything else
            from .algorithms import NextReactionMethod as Sampler

        return Sampler(self, state, tstart, tmax, steps)

    def sample(self, state, tstart=0., tmax=float('inf'), steps=None, seed=None):
        """Create trajectory sampler for given state

        The method returns an automatically chosen sampling algorithm
        suitable for the given stochastic process. The returned sampler
        can be iterated over to generate transitions along a trajectory:
 
        >>> process = Process()
        >>> state = {}
        >>> trajectory = process.trajectory(state)
        >>> for dt, transitions in trajectory():
        ...     print(trajectory.time, trajectory.state, transitions)

        C.f. stocal.algorithms.StochasticSimulationAlgorithm for details
        on the sampler class.
        """
        # This method yields transitions according to the future
        # StochasticSimulationAlgorithm specification in the form
        # (dt, transition_dct). It is planned to replace the current
        # trajectory method.
        sampler = self._trajectory(state, tstart=tstart, tmax=tmax, steps=steps, seed=seed)

        class _Wrapper(object):
            def __getattr__(self, attr):
                return getattr(sampler, attr)

            def __setattr__(self, attr, val):
                return setattr(sampler, attr, val)

            def __iter__(self):
                time = sampler.time
                for transition in sampler:
                    yield sampler.time-time, {transition: 1}
                    time = sampler.time
        return _Wrapper()

    def flatten(self, initial_species, max_steps=1000):
        Proc = type(self)
        flat_process = Proc(self.transitions)

        novel_species = multiset({s:float('inf') for s in initial_species })
        species = multiset()

        for step in range(max_steps):
            next_species = multiset()
            for rule in self.rules:
                for trans in rule.infer_transitions(novel_species, species):
                    flat_process.transitions.append(trans)
                    next_species += {s: float('inf') for s in trans.true_products if s not in species}
            species += novel_species
            if not next_species.domain:
                break
            else:
                novel_species = next_species
        else:
            raise ValueError("Flattening did not converge within %d steps" % max_steps)

        return flat_process


class ReactionRule(TransitionRule):
    """Deprecated. Identical to TransitionRule"""
    def __init__(self, *args, **opts):
        warnings.warn("Use TransitionRule instead", DeprecationWarning)
        super(ReactionRule, self).__init__(*args, **opts)

