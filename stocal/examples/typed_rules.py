"""An example on how to use typed rules
"""
import stocal

ABlock = type("ABlock", (str,), {})
BBlock = type("BBlock", (str,), {})


state = {ABlock(t): 100 for t in "xyz"}
state.update({BBlock(t): 100 for t in "uvw"})


class Polymerization(stocal.ReactionRule):
    """Alternate polymerization of ABlock's and BBlock's

    Polymerization only occurs among polymers of distinct type
    ABlock and BBlock. The type of the product polymer equals the
    right-most monomer type.
    """
    Transition = stocal.MassAction
    signature = [ABlock, BBlock]

    def novel_reactions(self, a, b):
        assert all(isinstance(species, typ) for species, typ in zip(
            (a, b), self.signature))
        yield self.Transition([a, b], [BBlock(a+b)], 1.)
        yield self.Transition([a, b], [ABlock(b+a)], 1.)


process = stocal.Process(rules=[Polymerization()])

if __name__ == '__main__':
    traj = process.trajectory(state, steps=1000)
    for trans in traj:
        print(traj.time, traj.state)
