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

"""A framework for rule-based stochastic simulation

This package offers functionalities to perform stochastic simulations
over regular as well as rule-based systems. To run a stochastic
simulation, simply define a Process as a collection of Transition and
Rule objects and call its trajectory method, providing an initial
state and a stop criterion (max simulated time or number of steps).
Iterating over the trajectory yields a series of fired transitions.
After ech transition, the attributes trajectory.time and
trajectory.state give information about the system state.

A user tutorial is provided in stocal/doc/tutorial.html
Example uses can be found in the stocal.examples package.
"""
import abc


class Transition(object) :
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
	__metaclass__ = abc.ABCMeta

	rule = None

	def __init__(self, reactants, products) :
		"""Initialization
		
		reactants and products are either mappings that give the
		stoichiometric factor of each involved species, or sequences
		which are interpreted as unordered lists.
		"""
		# if reactants/products are not a dictionary, they are assumed
		# to be unordered sequences of reactants/products
		def seq_to_dict(seq) :
			return seq if isinstance(seq, dict) else {
				s: seq.count(s) for s in set(seq)
			}
		reactants = seq_to_dict(reactants)
		products = seq_to_dict(products)

		if not all(n>0 for n in reactants.itervalues()) :
			raise ValueError("reactant stoichiometries must be positive.")
		if not all(n>0 for n in products.itervalues()) :
			raise ValueError("product stoichiometries must be positive.")
		if not reactants and not products :
			raise ValueError(
				"%s must have either reactants or products."
				% type(self).__name__
			)

		self.reactants = reactants
		self.products = products

		self.true_reactants = {
			s : n - products.get(s,0)
			for s,n in reactants.iteritems()
			if products.get(s,0) < n
		}
		self.true_products = {
			s : n - reactants.get(s,0)
			for s,n in products.iteritems()
			if reactants.get(s,0) < n
		}

		self.last_occurrence = -1.
		self._hash = 0

	def __eq__(self, other) :
		"""Structural congruence

		Transitions are equal if their reactants and products are equal.
		"""
		return (
			isinstance(other, Transition) and
			self.reactants == other.reactants and
			self.products == other.products
		)

	def __ne__(self, other) :
		return not self==other

	def __hash__(self) :
		if not self._hash :
			self._hash = hash((
				tuple(sorted(self.reactants.items())),
				tuple(sorted(self.products.items())),
			))
		return self._hash

	def __repr__(self) :
		try :
			return '%s(%s, %s)' % (
				type(self).__name__, self.reactants, self.products
			)
		except AttributeError :
			return super(Transition, self).__repr__()

	def __str__(self) :
		def dct2str(dct) :
			return ' + '.join(
				s if n==1 else '%s*%s' % (n,s) for s,n in dct.iteritems()
			)
		try :
			return '%s --> %s' % (dct2str(self.reactants), dct2str(self.products))
		except AttributeError :
			return super(Transition, self).__str__()

	@abc.abstractmethod
	def next_occurrence(self, time, state) :
		"""Time of next occurrence after given time.
		
		This method has to be implemented by a subclass. See
		Event.next_occurrence and Reaction.next_occurrence for examples.
		"""
		return float('inf')


class Reaction(Transition) :
	"""Abstract base class for stochastic transitions.

	Stochastic transitions are those that occur with a certain
	proponsity within a given time interval. 
	Subclasses must implement the propensity method.
	"""
	def __eq__(self, other) :
		"""Structural congruence
		
		Reactions are equal if their reactants, products, and propensity
		functions are equal.
		"""
		return (
			super(Reaction,self).__eq__(other) and
			isinstance(other, Reaction) and
			type(self).propensity == type(other).propensity
		)

	def __hash__(self) :
		if not self._hash :
			self._hash = hash((
				super(Reaction, self).__hash__(),
				type(self).propensity
			))
		return self._hash

	@abc.abstractmethod
	def propensity(self, state) :
		"""Reaction propensity.

		This method has to be provided by a subclass and must return a
		non-negative float denoting the reaction propensity for the
		provided state. The function should not modify the provided state.
		"""
		return 0.

	def next_occurrence(self, time, state) :
		"""Determine next reaction firing time.

		This is a helper function to use Reactions in next-firing-time
		based TrajectorySampler's. The method randomly draws a delay
		from a Poisson distribution with mean propensity  and returns
		the given current time plus the delay.
		"""
		from random import random
		from math import log
		p = self.propensity(state)
		if not p : return float('inf')
		else : return time - log(random())/p


class MassAction(Reaction) :
	"""Reactions with mass action kinetics.
	
	The propensity of a mass action reaction is defined as the
	stochastic rate constant c times the number of possible reaction
	partners, where the latter is the product of the binomial
	coefficients to choose n reaction partners out of m molecules for
	each reacting species.
	"""
	def __init__(self, reactants, products, c) :
		if c < 0 :
			raise ValueError("stochastic rate constants must be non-negative.")
		super(MassAction, self).__init__(reactants, products)
		self.c = c

	def __repr__(self) :
		try :
			return '%s(%s, %s, %g)' % (
				type(self).__name__, self.reactants, self.products, self.c
			)
		except AttributeError :
			return super(MassAction, self).__repr__()

	def propensity(self, state) :
		"""Reaction propensity for the given state.
		
		Calling propensity does not modify the underlying reaction.
		"""
		def choose(n,k) :
			return reduce(lambda x,i: x*(n+1-i)/i, xrange(1,k+1), 1)
		return reduce(
			lambda a,b: a*b,
			(choose(state.get(s,0), n) for s,n in self.reactants.iteritems()),
			self.c
		)

	def __eq__(self, other) :
		"""Structural congruence
		
		MassAction reactions are equal if their reactants, products,
		and stochastic rate constants are equal.
		"""
		return (
			super(MassAction,self).__eq__(other) and
			isinstance(other, MassAction) and
			self.c == other.c
		)

	def __hash__(self) :
		if not self._hash :
			self._hash = hash((
				super(MassAction, self).__hash__(),
				self.c
			))
		return self._hash


class Event(Transition) :
	"""Deterministic transitions.
	
	Events are Transition's that occur either once at a specified time,
	or periodically with a given frequency starting at a specified time.
	"""
	def __init__(self, reactants, products, time, frequency=0) :
		if time<0 :
			raise ValueError("time must be greater than 0.")
		if frequency<0 :
			raise ValueError("dt must be greater than (or equal to) 0.")
		super(Event, self).__init__(reactants, products)
		self.t = time
		self.dt = frequency

	def __eq__(self, other) :
		"""Structural congruence
		
		Events are equal if their reactants, products, time, and
		frequency are equal.
		"""
		return (
			super(Event,self).__eq__(other) and
			isinstance(other, Event) and
			self.t == other.t and
			self.dt == other.dt
		)

	def __hash__(self) :
		if not self._hash :
			self._hash = hash((
				super(Event, self).__hash__(),
				self.t, self.dt
			))
		return self._hash

	def next_occurrence(self, time, state={}) :
		"""Next occurrence of the Event at or after time.

		If the event event does not re-occur, returns float('inf').
		Calling next_occurrence leaves the event unmodified.
		"""
		if self.dt :
			t = time + (self.t-time)%self.dt
			return t if self.last_occurrence!=time else t+self.dt
		elif time < self.t :
			return self.t
		elif time == self.t and self.last_occurrence != time :
			return time
		else :
			return float('inf')


class Rule(object) :
	"""Abstract base class for rules
	
	Subclasses must provide a class attribute called Transition,
	which denotes the class of Transitions they generate.
	They also must provide a method infer_transitions which performs
	the actual transition inference.
	"""
	__metaclass__ = abc.ABCMeta
	
	@abc.abstractproperty
	def Transition(self) :
		"""The type of Transition that the rule generates"""
		return Transition

	@abc.abstractmethod
	def infer_transitions(self, new_species, state) :
		"""infer new transitions among new_species and state
		
		This method is called by the stochastic simulation algorithm
		"in the middle of a transition", i.e. after removing reactants
		from the state but before adding new_species as products.
		Implementations must return an iterable of Transition objects.
		"""
		# XXX Should this take last_transition?
		raise StopIteration


class ReactionRule(Rule) :
	"""Abstract base class that facilitates inference of Reactions

	This class provides a standard implementation of infer_transitions
	that generates all combinations of length self.order which containt
	spieces in state and or new_species, and which became possible only
	with the last transition, but where not possible using species of
	state alone. I.e. at least one molecule in the combination must
	come from new_species.
	The inference algorithm then calls ReactionRule.novel_transitions,
	passing in each novel combination as an unordered list. This method,
	to be implemented by a subclass, should return an iterable over
	every reaction that takes the novel species as reactants.
	"""

	@abc.abstractproperty
	def order(self) :
		"""Reaction order of infered reactions.
		
		The order of a reaction is the number of reactant molecules.
		To be defined by a subclass."""
		return int()

	@abc.abstractmethod
	def novel_reactions(self, *reactants) :
		"""Infer reactions for the given unordered list of reactants.

		To be implemented by a subclass.
		"""
		raise StopIteration

	def infer_transitions(self, last_products, state) :
		"""Standard inference algorithm for Reactions.
		
		see help(type(self)) for an explanation of the algorithm.
		"""
		def combinations(reactants, species, novel) :
			if len(reactants) == self.order :
				if novel :
					yield reactants
				return
			if not species : return
			s, n, end = species.pop(0)
			if not novel and n>end : return
			for combination in combinations(reactants, list(species), novel) :
				yield combination
			m = min(self.order-len(reactants), end)
			for i in xrange(1,m+1) :
				reactants.append(s)
				for combination in combinations(reactants, list(species), novel or i>=n) :
					yield combination
			del reactants[-m:]

		species = sorted((
			(s, state.get(s,0)+1, min(last_products.get(s,0)+state.get(s,0), self.order))
			for s in set(last_products).union(state)
		), key=lambda (s,n,m) : n-m)
		for reactants in combinations([], species, False) :
			for trans in self.novel_reactions(*reactants) :
				yield trans


class Process(object) :
	"""Stochastic process class

	A collection of all transitions and rules that define a
	stochastic process. When initializing a TrajectorySampler with
	a Process instance, transitions get copied over into the sampler.
	This makes it possible to use a single process instance with
	multiple samplers.
	"""
	def __init__(self, transitions=[], rules=[]) :
		self.transitions = transitions
		self.rules = rules

	def trajectory(self, state, t=0., tmax=float('inf'), steps=None) :
		"""Create trajectory sampler for given state
		
		If any static or infered transition is deterministic, this returns
		the FirstReactionMethod, otherwise the DirectMethod."""
		Sampler = DirectMethod if (
			all(isinstance(r, Reaction) for r in self.transitions)
			and all(issubclass(r.Transition, Reaction) for r in self.rules)
		) else FirstReactionMethod

		return Sampler(self, state, t, tmax, steps)


class TrajectorySampler(object) :
	"""Abstract base class for stochastic trajectory samplers.

	This is the interface for stochastic trajectory samplers, i.e.
	imlementations of the stochastic simulation algorithm.
	A trajectory sampler is initialized with a given process and
	initial state, and options start time, end time, and maximal
	number of iterations. The sampler instance can then be iterated
	over to produce a stochastic trajectory:
	
	>>> trajectory = DirectMethod(process, state, steps=10000)
	>>> for transition in trajectory :
	...     print trajectory.time, trajectory.state, transition
	"""
	__metaclass__ = abc.ABCMeta
	
	def __init__(self, process, state, t=0., tmax=float('inf'), steps=None) :
		"""Initialize the sampler for the given Process and state.
		
		State is a dictionary that maps chemical species to positive
		integers denoting their copy number. An optional start time
		t, end time tmax, and maximal number of steps can be provided.
		"""
		if t<0 :
			raise ValueError("t must not be negative.")
		if tmax<0 :
			raise ValueError("tmax must not be negative.")
		if steps is not None and steps<0 :
			raise ValueError("steps must not be negative.")
		if any(n<=0 for n in state.itervalues()) :
			raise ValueError("state copy numbers must be positive")

		self.process = process
		self.step = 0
		self.steps = steps
		self.time = t
		self.tmax = tmax

		self.transitions = []
		for transition in self.process.transitions :
			self.add_transition(transition)
		self.state = {}
		self.update_state(state)

	@abc.abstractmethod
	def add_transition(self, transition) :
		"""Add a new transition to the sampler"""
		pass

	@abc.abstractmethod
	def update_state(self, dct, **kwds) :
		"""Modify sampler state.
		
		New transitions get infered appropriately."""
		pass

	@abc.abstractmethod
	def __iter__(self) :
		"""Sample stochastic trajectory.
		
		This method iteratively picks an applicable transition with
		probability proportional to its propensity, updates the sampler
		state accordingly, and yields the transition to the caller.
		After each transition, rules are invoked to infer potential
		novel transitions.		
		"""
		raise StopIteration
		yield None

	def _perform_transition(self, transition) :
		"""Perform the given transition.
		
		Remove reactants from and add products to state, and
		call self.rules to potentially infer novel transitions
		for the changed system state.
		"""
		def begin() :
			for species,n in transition.true_reactants.iteritems() :
				self.state[species] -= n
				if not self.state[species] :
					del self.state[species]
		def end() :
			for species,n in transition.true_products.iteritems() :
				if species not in self.state :
					self.state[species] = 0
				self.state[species] += n
		begin()
		for rule in self.process.rules :
			# XXX determine from state which rules should be queried
			for trans in rule.infer_transitions(transition.true_products, self.state) :
				trans.rule = rule
				self.add_transition(trans)
		end()


class DirectMethod(TrajectorySampler) :
	"""Implementation of Gillespie's direct method.
	
	The standard stochastic simulation algorithm, published in
	D. T. Gillespie, J. Comp. Phys. 22, 403-434 (1976).

	DirectMethod works only over processes that employ stochastic
	transitions (Reaction's). If the process includes Event's,
	FirstReactionMethod has to be used instead.

	See help(TrajectorySampler) for usage information.
	"""
	def __init__(self, process, state, t=0., tmax=float('inf'), steps=None) :
		super(DirectMethod, self).__init__(process, state, t, tmax, steps)
		if any(not isinstance(r,Reaction) for r in process.transitions) :
			raise ValueError("DirectMethod only works with Reactions.")
		if any(not issubclass(r.Transition, Reaction) for r in process.rules) :
			raise ValueError("DirectMethod only works with Reactions.")

	def add_transition(self, transition) :
		self.transitions.append(transition)

	def update_state(self, dct) :
		for rule in self.process.rules :
			for r in rule.infer_transitions(dct, self.state) :
				if r not in self.transitions :
					r.rule = rule
					self.add_transition(r)
		self.state.update(dct)
		for species in dct :
			if not species :
				del self.state[species]

	def __iter__(self) :
		from random import random
		from math import log
		from itertools import izip

		while True :
			if self.steps is not None and self.step==self.steps : return
			if self.time>=self.tmax : break

			propensities = [r.propensity(self.state) for r in self.transitions]

			total_propensity = sum(propensities)
			if not total_propensity : break

			# housekeeping: remove depleted transitions
			depleted = [
				i for i,(p,t) in enumerate(izip(propensities,self.transitions))
				if p==0. and t.rule
			]
			for i in reversed(depleted) :
				del self.transitions[i]
				del propensities[i]

			dt = -log(random())/total_propensity
			self.time += dt

			if self.time >= self.tmax :
				break
			else :
				pick = random()*total_propensity
				for p,transition in izip(propensities, self.transitions) :
					pick -= p
					if pick < 0. : break
				self.step += 1
				self._perform_transition(transition)
				yield transition

		if self.tmax<float('inf') : self.time = self.tmax


class FirstReactionMethod(TrajectorySampler) :
	"""Implementation of Gillespie's first reaction method.
	
	A stochastic simulation algorithm, published in
	D. T. Gillespie, J. Comp. Phys. 22, 403-434 (1976).

	FirstReactionMethod works with processes that feature deterministic
	transitions, i.e. Event's.

	See help(TrajectorySampler) for usage information.
	"""
	def __iter__(self) :
		while True :
			if self.steps is not None and self.step==self.steps : return
			if self.time>=self.tmax : break

			firings = [
				(trans.next_occurrence(self.time, self.state), trans)
				for trans in self.transitions
			]

			if not firings : break
			time, transition = min(firings)

			# housekeeping: remove depleted transitions
			depleted = [
				i for i,(t,r) in enumerate(firings)
				if t==float('inf') and (r.rule or isinstance(r, Event))
			]
			for i in reversed(depleted) :
				del self.transitions[i]
				del firings[i]

			if time >= self.tmax :
				break
			else :
				self.step += 1
				self.time = time
				transition.last_occurrence = time
				self._perform_transition(transition)
				yield transition

		if self.tmax<float('inf') : self.time = self.tmax

	def add_transition(self, transition) :
		self.transitions.append(transition)

	def update_state(self, dct) :
		for rule in self.process.rules :
			for r in rule.infer_transitions(dct, self.state) :
				if r not in self.transitions :
					r.rule = rule
					self.add_transition(r)
		self.state.update(dct)
		for species in dct :
			if not species :
				del self.state[species]
