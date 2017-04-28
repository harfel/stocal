import abc

class Reaction(object) :
	"""Reaction class

	Base class for stochastic transitions of reactants into products.
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
		
		This method has to be provided by a subclass.
		"""
		return 0.


class MassAction(Reaction) :
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
