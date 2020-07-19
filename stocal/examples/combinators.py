"""Stochastic combinator logic 

This example implements an aritficial chemistry of combinator logic expressions. Combinator logic
is equivalent to the lambda calculus but does not require variables. Expressions are binary trees
of the following combinators (where left associativity, CDE=(CD)E is implied):
    the identify combinator:     Ix... --> x...
    the constant combinator:     Kxy... --> x...
    the substitution combinator: Sfgx... --> fx(gx)...

Kruszewski and Mikolov introduced a Combinarory Chemistry in Proc. ALIFE Conf 2020
(https://doi.org/10.1162/isal_a_00258) as a combinator-conservative variant of the above:
          Ix... --> x... + I
         Kxy... --> x... + K + y
    x + Sgfx... --> fx(gx)... + S
"""
import stocal

I = 'I'
S = 'S'
K = 'K'
n = 1

def norm(first, *rest):
    while isinstance(first, tuple):
        rest = first[1:]+rest
        first = first[0]
    return (first,) + rest if rest else first


class Join(stocal.TransitionRule):
    Transition = stocal.MassAction
    c = 0.0001/n

    def novel_reactions(self, k, l):
        yield self.Transition([k, l], [norm(k,l)], self.c)
        yield self.Transition([k, l], [norm(l,k)], self.c)


class CombinatorRule(stocal.TransitionRule):
    Transition = stocal.MassAction
    c = 1.

    def all_redexes(self, term, react=None):
        try:
            yield self.react(react, *term)
        except TypeError:
            pass

        for i,x in enumerate(term):
            if isinstance(x, tuple):
                for y, prod in self.all_redexes(x):
                    yield term[:i] + y + term[i+1:], prod

    def novel_reactions(self, term):
        for y, prod in self.all_redexes(term):
            yield self.Transition([term], [norm(y)]+prod, self.c)


class ReduceI(CombinatorRule):
    def react(self, react, combinator, x, *rest):
        if combinator != I:
            raise TypeError
        return (x,) + rest, [I]

class ReduceK(CombinatorRule):
    def react(self, react, combinator, x, y, *rest):
        if combinator != K:
            raise TypeError
        return (x,) + rest, [K, y]

class ReduceS(CombinatorRule):
    c = 1.*CombinatorRule.c/n

    def react(self, react, combinator, f, g, x, *rest):
        if combinator != S:
            raise TypeError
        if x != react:
            raise TypeError
        return (f,x,(g,x))+rest, [S]

    def novel_reactions(self, left, right):
        for y, prod in self.all_redexes(left, right):
            yield self.Transition([left, right], [norm(y)]+prod, self.c)
        for y, prod in self.all_redexes(right, left):
            yield self.Transition([left, right], [norm(y)]+prod, self.c)


process = stocal.Process(rules=[ReduceI(), ReduceK(), ReduceS(), Join()])
initial_state = {I: n, K: n, S: n}

if __name__ == '__main__':
    traj = process.sample(initial_state, steps=100000)
    for t,trans in traj:
        print(trans.keys()[0].rule.__class__.__name__, str(trans.keys()[0]),traj.time, traj.state)
    print(traj.state)
