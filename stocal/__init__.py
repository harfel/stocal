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


class Reaction(Transition) :
	"""Reaction class

	Reactions are stochastic Transitions.
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


class MassAction(Reaction) :
	"""Reaction with mass action kinetics"""
	def __init__(self, reactants, products, c) :
		if c < 0 :
			raise ValueError("stochastic rate constants must be non-negative.")
		super(MassAction, self).__init__(reactants, products)
		self.c = c

	def propensity(self, state) :
		def choose(n,k) :
			return reduce(lambda x,i: x*(n+1-i)/i, xrange(1,k+1), 1)
		return reduce(
			lambda a,b: a*b,
			(choose(state.get(s,0), n) for s,n in self.reactants.iteritems()),
			self.c
		)


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

	def trajectory(self, state, t=0., tmax=-1., steps=-1) :
		"""Create trajcetory sampler for given state"""
		return DirectMethod(self, state, t, tmax, steps)


class TrajectorySampler(object) :
	__metaclass__ = abc.ABCMeta
	
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


from random import random
from math import log
from itertools import izip

class DirectMethod(TrajectorySampler) :
	"""Implementation of Gillespie's direct method"""
	def __init__(self, process, state, t=0., tmax=-1., steps=-1) :
		if t<0 :
			raise ValueError("t must not be negative.")
		if any(not isinstance(r,Reaction) for r in process.transitions) :
			raise ValueError("DirectMethod only works with Reactions.")
		if any(not issubclass(r.Transition, Reaction) for r in process.rules) :
			raise ValueError("DirectMethod only works with Reactions.")
		self.process = process
		self.state = state
		self.transitions = []
		self.step = 0
		self.steps = steps
		self.time = t
		self.tmax = tmax

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
		while True :
			if self.steps>0 and self.step==self.steps : break
			if self.tmax>=0 and self.time>=self.tmax : break

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
				if self.tmax>=0 : self.time = self.tmax
				break

			dt = -log(random()/total_propensity)
			self.time += dt
			self.step += 1

			if self.time >= self.tmax >= 0 :
				self.time = self.tmax
				break
			else :
				pick = random()*total_propensity
				for p,transition in izip(propensities, self.transitions) :
					pick -= p
					if pick < 0. : break
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
