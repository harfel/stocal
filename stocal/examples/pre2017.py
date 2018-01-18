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
	Transition = stocal.MassAction
	c = 1.

	def __str__(self) :
		return 'kl --> k+l'

	def infer_transitions(self, new_species, state) :
		for kl in new_species :
			if kl in state : continue
			for i in xrange(1,len(kl)) :
				k = kl[:i]
				l = kl[i:]
				f = 2 if kl[-i:]==k and kl[:-i]==l and 2*i!=len(kl) else 1
				yield self.Transition([kl], [k, l], f*self.c)


class LigationRule(stocal.Rule) :
	Transition = stocal.MassAction
	c = 1.e-10

	def __str__(self) :
		return 'k+l --> kl'

	def infer_transitions(self, new_species, state) :
		for k in new_species :
			if state.get(k,0) >= 2 : continue
			for l in new_species :
				f = 2 if k+l==l+k else 1
				yield self.Transition([k,l], [k+l], f*self.c)
			for l in state :
				f = 2 if k+l==l+k else 1
				yield self.Transition([k,l], [k+l], f*self.c)
				yield self.Transition([k,l], [l+k], f*self.c)


class AutoCatalysisRule(stocal.Rule) :
	Transition = stocal.MassAction
	c = 1.e-6

	def __str__(self) :
		return 'k+l+kl --> 2*kl'

	def infer_transitions(self, new_species, state) :
		for k in new_species :
			if state.get(k,0) >= 2 : continue
			for l in new_species :
				if k+l in state :
					f = 2 if k+l==l+k else 1
					yield self.Transition([k,l,k+l], {k+l:2}, f*self.c)
				if k.startswith(l) and k!=l and l in state :
					f = 2 if k[len(l):]+l==l+k[len(l):] else 1
					yield self.Transition([k[len(l):],l,k], {k:2}, f*self.c)
			for l in state :
				if k+l in state :
					f = 2 if k+l==l+k else 1
					yield self.Transition([k,l,k+l], {k+l:2}, f*self.c)
				if l+k in state :
					f = 2 if k+l==l+k else 1
					yield self.Transition([k,l,l+k], {l+k:2}, f*self.c)
				if k.startswith(l) and k!=l :
					f = 2 if k[len(l):]+l==l+k[len(l):] else 1
					yield self.Transition([k[len(l):],l,k], {k:2}, f*self.c)
				if k.endswith(l) and k!=l :
					f = 2 if k[:-len(l)]+l==l+k[:-len(l)] else 1
					yield self.Transition([k[:-len(l)],l,k], {k:2}, f*self.c)


process = stocal.Process(rules=[DegradationRule(), LigationRule(), AutoCatalysisRule()])
initial_state = { 'a': 200000, 'b': 200000 }


if __name__ == '__main__' :
	traj = stocal.DirectMethod(process, initial_state, tmax=1)
	for _ in traj :
		print(traj.time, traj.state)
