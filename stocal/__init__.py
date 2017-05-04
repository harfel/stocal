import abc

class Transition(object) :
	"""Transitions are transformations of reactants into products

	After initialization, reactants and products must not be altered.
	Doing so leads to undefined behavior.
	"""
	__metaclass__ = abc.ABCMeta

	rule = None

	def __init__(self, reactants, products) :
		"""Initialization
		
		reactants and products are mappings that give the stoichiometric
		factor of each involved species.
		"""
		if not all(n>0 for n in reactants.itervalues()) :
			raise ValueError("reactant stoichiometries must be positive.")
		if not all(n>0 for n in products.itervalues()) :
			raise ValueError("product stoichiometries must be positive.")
		if not reactants and not products :
			raise ValueError(
				"%s must have either reactants or products." % type(self).__name__
			)

		self.reactants = reactants
		self.products = products
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

	def __hash__(self) :
		if not self._hash :
			self._hash = hash((
				tuple(sorted(self.reactants.items())),
				tuple(sorted(self.products.items())),
			))
		return self._hash

	def __repr__(self) :
		return '%s(%s, %s)' % (
			type(self).__name__, self.reactants, self.products
		)

	def __str__(self) :
		def dct2str(dct) :
			return ' + '.join(
				s if n==1 else '%s*%s' % (n,s) for s,n in dct.iteritems()
			)
		return '%s --> %s' % (dct2str(self.reactants), dct2str(self.products))


class Reaction(Transition) :
	"""Stochastic Transitions

	Subclasses must implement a propensity method.
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
		"""Reaction propensity
		
		This method has to be provided by a subclass and must return a
		non-negative float denoting the reaction propensity for the
		provided state. The function should not modify the provided state.
		"""
		return 0.

	def next_ocurrence(self, time, state) :
		"""Determine next reaction firng
		
		This is a helper function to use Reactions in place of Events.
		It randomly draws a daly from a Poisson distribution and
		returns the given current time plus the delay.
		"""
		from random import random
		from math import log
		p = self.propensity(state)
		if not p : return float('inf')
		else : return time - log(random()/p)


class MassAction(Reaction) :
	"""Reactions with mass action kinetics"""
	def __init__(self, reactants, products, c) :
		if c < 0 :
			raise ValueError("stochastic rate constants must be non-negative.")
		super(MassAction, self).__init__(reactants, products)
		self.c = c

	def __repr__(self) :
		return '%s(%s, %s, %f)' % (
			type(self).__name__, self.reactants, self.products, self.c
		)

	def propensity(self, state) :
		def choose(n,k) :
			return reduce(lambda x,i: x*(n+1-i)/i, xrange(1,k+1), 1)
		return reduce(
			lambda a,b: a*b,
			(choose(state.get(s,0), n) for s,n in self.reactants.iteritems()),
			self.c
		)


class Event(Transition) :
	"""Deterministic transitions.
	
	Events are transitions that occur either once, at a specified time, or
	periodically with a given frequency starting at a specified time."""
	def __init__(self, reactants, products, time, frequency=0) :
		if frequency<0 :
			raise ValueError("dt must be greater than (or equal to) 0.")
		super(Event, self).__init__(reactants, products)
		self.t = time
		self.dt = frequency

	def next_occurrence(self, time, state={}) :
		if self.dt :
			return self.t + self.dt*((time-self.t)//self.dt+1)
		elif time < self.t :
			return self.t
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
		
		The method must return an iterable of Transition objects.
		"""
		raise StopIteration
		yield None


class Process(object) :
	"""Stochastic process class

	A collection of all transitions and rules that define a
	stochastic process.
	"""
	def __init__(self, transitions=[], rules=[]) :
		self.transitions = transitions
		self.rules = rules

	def trajectory(self, state, t=0., tmax=float('inf'), steps=-1) :
		"""Create trajcetory sampler for given state
		
		If any static or infered transition is deterministic, this returns
		the FirstReactionMethod, otherwise the DirectMethod."""
		Sampler = DirectMethod if (
			all(isinstance(r, Reaction) for r in self.transitions)
			and all(issubclass(r.Transition, Reaction) for r in self.rules)
		) else FirstReactionMethod

		return Sampler(self, state, t, tmax, steps)


class TrajectorySampler(object) :
	__metaclass__ = abc.ABCMeta
	
	def __init__(self, process, state, t=0., tmax=float('inf'), steps=-1) :
		if t<0 :
			raise ValueError("t must not be negative.")
		self.process = process
		self.step = 0
		self.steps = steps
		self.time = t
		self.tmax = tmax

	@abc.abstractmethod
	def add_transition(self, transition) :
		"""Add a new transition to the sampler"""
		pass

	@abc.abstractmethod
	def update_state(self, dct, **kwds) :
		"""Modify sampler state"""
		pass

	@abc.abstractmethod
	def __iter__(self) :
		"""Iteratively apply and return a firing transition"""
		raise StopIteration
		yield None


class DirectMethod(TrajectorySampler) :
	"""Implementation of Gillespie's direct method"""
	def __init__(self, process, state, t=0., tmax=float('inf'), steps=-1) :
		super(DirectMethod, self).__init__(process, state, t, tmax, steps)
		if any(not isinstance(r,Reaction) for r in process.transitions) :
			raise ValueError("DirectMethod only works with Reactions.")
		if any(not issubclass(r.Transition, Reaction) for r in process.rules) :
			raise ValueError("DirectMethod only works with Reactions.")
		self.state = state
		self.transitions = []
		for transition in process.transitions :
			self.add_transition(transition)

	def add_transition(self, transition) :
		"""Add a new transition to the sampler"""
		self.transitions.append(transition)

	def update_state(self, dct) :
		"""Modify sampler state"""
		for rule in self.process.rules :
			for r in rule.infer_transitions(dct, self.state) :
				if r not in self.transitions :
					r.rule = rule
					self.add_transition(r)
		self.state.update(dct)

	def __iter__(self) :
		"""Iteratively apply and return a firing transition"""
		from random import random
		from math import log
		from itertools import izip

		while True :
			if self.steps>0 and self.step==self.steps : break
			if self.time>=self.tmax : break

			propensities = [r.propensity(self.state) for r in self.transitions]

			# housekeeping: remove depleted transitions
			depleted = [
				i for i,p in enumerate(izip(propensities,self.transitions))
				if p==0. and r.rule
			]
			for i in reversed(depleted) :
				del self.transitions[i]
				del propensities[i]

			total_propensity = sum(propensities)
			if not total_propensity :
				if self.tmax<float('inf') : self.time = self.tmax
				break

			dt = -log(random()/total_propensity)
			self.time += dt

			if self.time >= self.tmax :
				self.time = self.tmax
				break
			else :
				pick = random()*total_propensity
				for p,transition in izip(propensities, self.transitions) :
					pick -= p
					if pick < 0. : break
				self.step += 1
				self._perform_transition(transition)
				yield transition

	def _perform_transition(self, transition) :
		def begin() :
			for species,n in transition.reactants.iteritems() :
				self.state[species] -= n
				if not self.state[species] :
					del self.state[species]
		def end() :
			for species,n in transition.products.iteritems() :
				if species not in self.state :
					self.state[species] = 0
				self.state[species] += n
		begin()
		for rule in self.process.rules :
			for r in rule.infer_transitions(transition.products, self.state) :
				if r not in self.transitions :
					r.rule = rule
					self.add_transition(r)
		end()


class FirstReactionMethod(TrajectorySampler) :
	def __init__(self, process, state, t=0., tmax=float('inf'), steps=-1) :
		super(FirstReactionMethod, self).__init__(process, state, t, tmax, steps)
		self.state = state
		self.transitions = []
		for transition in process.transitions :
			self.add_transition(transition)

	def __iter__(self) :
		"""Iteratively apply and return a firing transition"""
		while True :
			if self.steps>0 and self.step==self.steps : break
			if self.time>=self.tmax : break

			firings = [
				(r.next_ocurrence(self.time, self.state),r)
				for r in self.transitions
			]

			# housekeeping: remove depleted transitions
			depleted = [
				i for i,(p,r) in enumerate(firings)
				if p==float('inf') and (r.rule or isinstance(r, Event))
			]
			for i in reversed(depleted) :
				del self.transitions[i]
				del firings[i]

			if not firings : break

			time, transition = min(firings)

			if time >= self.tmax :
				break
			else :
				self.step += 1
				self.time = time
				self._perform_transition(transition)
				yield transition

		if self.tmax<float('inf') : self.time = self.tmax

	def add_transition(self, transition) :
		"""Add a new transition to the sampler"""
		self.transitions.append(transition)

	def update_state(self, dct) :
		"""Modify sampler state"""
		for rule in self.process.rules :
			for r in rule.infer_transitions(dct, self.state) :
				if r not in self.transitions :
					r.rule = rule
					self.add_transition(r)
		self.state.update(dct)

	def _perform_transition(self, transition) :
		def begin() :
			for species,n in transition.reactants.iteritems() :
				self.state[species] -= n
				if not self.state[species] :
					del self.state[species]
		def end() :
			for species,n in transition.products.iteritems() :
				if species not in self.state :
					self.state[species] = 0
				self.state[species] += n
		begin()
		for rule in self.process.rules :
			for r in rule.infer_transitions(transition.products, self.state) :
				if r not in self.transitions :
					r.rule = rule
					self.add_transition(r)
		end()
