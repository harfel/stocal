"""An example on how to use typed rules
"""
import stocal

ABlock = stocal.molecular_type("ABlock")
BBlock = stocal.molecular_type("ABlock")


state = {
    ABlock('x'): 100,
    ABlock('y'): 100,
    ABlock('z'): 100,
    BBlock('u'): 100,
    BBlock('v'): 100,
    BBlock('w'): 100,
}


class Polymerization(stocal.ReactionRule):
    """Alternate polymerization of ABlock's and BBlock's

    Polymerization only occurs among polymers of distinct type
    ABlock and BBlock. The type of the product polymer equals the
    right-most monomer type.
    """
    Transition = stocal.MassAction
    signature = [ABlock, BBlock]

    def novel_reactions(self, a, b):
        yield self.Transition([a, b], [BBlock(a+b)], 1.)
        yield self.Transition([a, b], [ABlock(b+a)], 1.)


process = stocal.Process(rules=[Polymerization()])

if __name__ == '__main__':
    traj = process.trajectory(state, steps=1000)
    for trans in traj:
        print(traj.time, traj.state)
