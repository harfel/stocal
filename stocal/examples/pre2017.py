"""Exact autocatalytic polymer model

This model implements the system published in

    H Fellermann, S Tanaka, R Rasmussen, Sequence selection by
    dynamical symmetry breaking in an autocatalytic binary polymer
    model, Phys. Rev. E 96, 062407 (2017)
    https://doi.org/10.1103/PhysRevE.96.062407

The model is based on three rules that determine the behavior of
self-replicating heteropolymers:

  * The DegradationRule allows a polymer to break into two
    substrands (k.l --> k + l)

  * The LigationRule allows two random polymers to form a
    non-commutative ligation product (k+ l --> k.l)

  * The AutoCatalysisRule allows a polymer to be replicated
    catalytically from two exactly matching substrands
    (k + l + k.l --> 2*k.l)

The rules operate over any monomer alphabet, but the model defines
an initial state with two types of monomers, making this a binary exact
autocatalytic polymer model.
"""

import stocal


class DegradationRule(stocal.Rule) :
	class Transition(stocal.MassAction) :
		c = 1.

		def __init__(self, kl, i) :
			super(DegradationRule.Transition, self).__init__(
				{kl: 1},
				{kl[:i]: 2} if kl[:i]==kl[i:] else {kl[:i]: 1, kl[i:]: 1},
				2*self.c if kl[:i]==kl[i:] else self.c
			)

	def __str__(self) :
		return 'kl --> k+l'

	def infer_transitions(self, new_species, state) :
		for kl in new_species :
			if kl in state : continue
			for i in xrange(1,len(kl)) :
				yield self.Transition(kl, i)


class LigationRule(stocal.Rule) :
	class Transition(stocal.MassAction) :
		c = 1.e-10

		def __init__(self, k, l) :
			super(LigationRule.Transition, self).__init__(
				{k: 2} if k==l else {k: 1, l: 1},
				{k+l: 1},
				2*self.c if k==l else self.c
			)

	def __str__(self) :
		return 'k+l --> kl'

	def infer_transitions(self, new_species, state) :
		for k in new_species :
			if state.get(k,0) >= 2 : continue
			for l in new_species :
				yield self.Transition(k,l)
			for l in state :
				yield self.Transition(k,l)
				yield self.Transition(l,k)


class AutoCatalysisRule(stocal.Rule) :
	class Transition(stocal.MassAction) :
		c = 1.e-6

		def __init__(self, k, l) :
			super(AutoCatalysisRule.Transition, self).__init__(
				{k: 2, k+l: 1} if k==l else {k: 1, l: 1, k+l: 1},
				{k+l: 2},
				2*self.c if k==l else self.c
			)

	def __str__(self) :
		return 'k+l+kl --> 2*kl'

	def infer_transitions(self, new_species, state) :
		for k in new_species :
			if state.get(k,0) >= 2 : continue
			for l in new_species :
				if k+l in state :
					yield self.Transition(k,l)
				if k.startswith(l) and l in state :
					yield self.Transition(l, k[len(l):])
			for l in state :
				if k+l in state :
					yield self.Transition(k,l)
				if l+k in state :
					yield self.Transition(l,k)
				if k.startswith(l) and l in state :
					yield self.Transition(l, k[len(l):])
				if k.endswith(l) and l in state :
					yield self.Transition(k[:len(l)],l)


process = stocal.Process(rules=[DegradationRule(), LigationRule(), AutoCatalysisRule()])
initial_state = { 'a': 200000, 'b': 200000 }
