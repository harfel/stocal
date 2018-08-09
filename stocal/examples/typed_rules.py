"""An example of how to use typed rules

We model polymers made up of two types of monomers, with the contraint
that polymerization can only occurr among differently typed monomers.
We model this with a type system where the type of a polymer expresses
its leftmost and rightmost monomer types. E.g. a polymer starting with
A and ending in B has the type AB. (Monomers are therefor polymers of
type AA or BB.)

This example stretches a bit the boundaries of what can be done in a 
type system.
"""
import stocal

# We first define types for the polymers that can occurr in our chemistry

AA = stocal.molecular_type("AA")
AB = stocal.molecular_type("AB")
BA = stocal.molecular_type("BA")
BB = stocal.molecular_type("BB")

# We now define a base class for generic polymerzation. We will derive
# typed polymerization rules from this genereic base class.

class Polymerization(stocal.TransitionRule):
    """Generic polymerization rule.
    
    Polymerization.Transition's are MassAction reactions.
    A polymerization of two reactants can form either one or two
    products depending on the types of the leftmost and rightmost
    monomers. For example, AA and BB can form AB or BA. These two
    product types are stored in the attributes Polymeriation.KL
    and Polymerization.LK which are to be overloaded by subclasses.
    If either of them is None (default), reactants are not allowed
    to bind in this way.
    """
    Transition = stocal.MassAction
    KL = None
    LK = None

    def novel_reactions(self, k, l):
        """Infer reactions with currect product types."""
        if self.KL:
            yield self.Transition([k, l], [self.KL(k+l)], 1.)
        if self.LK:
            yield self.Transition([l, k], [self.LK(l+k)], 1.)


# We now define typed specializations of all possible Polymerization
# rules in our chemistry. Each rule has a typed signature for its
# reactants and defines one or two product types.

class AA_AB(Polymerization):
    signature = [AA, AB]
    LK = AA

class AA_BA(Polymerization):
    signature = [AA, BA]
    KL = AA

class AA_BB(Polymerization):
    signature = [AA, BB]
    KL = AB
    LK = BA

class AB_AB(Polymerization):
    signature = [AB, AB]
    KL = AB
    LK = AB

class AB_BB(Polymerization):
    signature = [AB, BB]
    LK = BB

class BA_BA(Polymerization):
    signature = [BA, BA]
    KL = BA
    LK = BA

class BA_BB(Polymerization):
    signature = [BA, BB]
    KL = BB

# The initial condition is just a bunch of monomers of types AA and BB

state = {
    AA('a'): 100,
    AA('b'): 100,
    AA('c'): 100,
    BB('x'): 100,
    BB('y'): 100,
    BB('z'): 100,
}

# The process is defined to use all subclasses of Polymerization

process = stocal.Process(rules=[R() for R in Polymerization.__subclasses__()])

# And we go in to sample the trajectory of the process...

if __name__ == '__main__':
    traj = process.sample(state, steps=1000)
    for trans in traj:
        print(traj.time, traj.state)
