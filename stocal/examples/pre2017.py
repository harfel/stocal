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

With the default parameters, simulations will take a significant
amount of time.
"""
import stocal


alpha = 1.e-10
beta = 1000**-2
initial_state = {c: 200000 for c in 'ab'}


class DegradationRule(stocal.ReactionRule):
    """Break a string into any two nonempty substrings"""
    Transition = stocal.MassAction

    def __str__(self):
        return 'kl --> k+l'

    def novel_reactions(self, kl):
        for i in range(1, len(kl)):
            k = kl[:i]
            l = kl[i:]
            yield self.Transition([kl], [k, l], 1.)


class LigationRule(stocal.ReactionRule):
    """Join any two strings into their concatenations"""
    Transition = stocal.MassAction

    def __str__(self):
        return 'k+l --> kl'

    def novel_reactions(self, k, l):
        yield self.Transition([k, l], [k+l], alpha)
        yield self.Transition([k, l], [l+k], alpha)


class AutoCatalysisRule(stocal.ReactionRule):
    """Replicate any string from two matching substrings"""
    Transition = stocal.MassAction

    def __str__(self):
        return 'k+l+kl --> 2*kl'

    def novel_reactions(self, k, l, m):
        k, l, m = sorted([k, l, m], key=len)
        if k+l == m:
            yield self.Transition([k, l, m], {m: 2}, beta)
        if l+k == m:
            yield self.Transition([k, l, m], {m: 2}, beta)


process = stocal.Process(
    rules=[DegradationRule(), LigationRule(), AutoCatalysisRule()]
)


if __name__ == '__main__':
    traj = process.trajectory(initial_state, tmax=100.)
    for _ in traj:
        print(traj.time, traj.state)
