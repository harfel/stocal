import abc

class Reaction(object) :
	"""Reaction class

	Abstract base class for stochastic transitions of reactants into products.
	Subclasses must implement the propensity method.

	After initialization, reactants and products must not be altered.
	Doing so leads to undefined behavior.
	"""
	__metaclass__ = abc.ABCMeta

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
			raise ValueError("Reaction must have either reactants or products.")

		self.reactants = reactants
		self.products = products
		self._hash = 0

	def __eq__(self, other) :
		"""Structural congruence
		
		Reactions are equal if their reactants, products, and propensity
		functions are equal.
		"""
		return (
			self.reactants == other.reactants and
			self.products == other.products and
			type(self).propensity == type(other).propensity
		)

	def __hash__(self) :
		if not self._hash :
			self._hash = hash((
				tuple(sorted(self.reactants.items())),
				tuple(sorted(self.products.items())),
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
	
	Rule subclasses must provide the infer_reactions method.
	"""
	__metaclass__ = abc.ABCMeta
	
	@abc.abstractmethod
	def infer_reactions(self, new_species, state) :
		raise StopIteration


class Process(object) :
	"""Stochastic process class

	A collection of all reactions (and later rules) that define a
	stochastic process.
	"""
	def __init__(self, reactions) :
		self.reactions = reactions

	def trajectory(self, state, t=0., tmax=-1., steps=-1) :
		"""Create trajcetory sampler for given state"""
		return DirectMethod(self, state, t, tmax, steps)


class TrajectorySampler(object) :
	__metaclass__ = abc.ABCMeta
	
	@abc.abstractmethod
	def add_reaction(self, reaction) :
		"""Add a new reaction to the sampler"""
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
		self.state = state
		self.reactions = []
		self.step = 0
		self.steps = steps
		self.time = t
		self.tmax = tmax

		for reaction in process.reactions :
			self.add_reaction(reaction)

	def add_reaction(self, reaction) :
		"""Add a new reaction to the sampler"""
		self.reactions.append(reaction)

	def update_state(self, dct, **kwds) :
		"""Modify sampler state"""
		self.state.update(dct, **kwds)

	def __iter__(self) :
		"""Iteratively apply and return a firing transition"""
		while True :
			if self.steps>0 and self.step==self.steps : break
			if self.tmax>=0 and self.time>=self.tmax : break
			propensities = [r.propensity(self.state) for r in self.reactions]
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
				for p,reaction in izip(propensities, self.reactions) :
					pick -= p
					if pick < 0. : break
				self._perform_reaction(reaction)
				yield reaction

	def _perform_reaction(self, reaction) :
		for species,n in reaction.reactants.iteritems() :
			self.state[species] -= n
			if not self.state[species] :
				del self.state[species]
		for species,n in reaction.products.iteritems() :
			if species not in self.state :
				self.state[species] = 0
			self.state[species] += n
