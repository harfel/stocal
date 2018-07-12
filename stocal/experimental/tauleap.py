"""tau leaping

This module provides an approximate stochastic simulation algorithm
based on tau-leaping, as derived in

Efficient step size selection for the tau-leaping simulation method
Y. Cao, D. T. Gillespie, L. R. Petzold, J. Chem. Phys. 124, 044109 (2006)

"""
import sys
import logging
from math import log
try:
    from numpy.random import RandomState
except ImportError:
    logging.error("stocal.experimental.tauleap requires numpy.")
    sys.exit(1)

from stocal.algorithms import DirectMethod
from stocal.structures import multiset


class CaoMethod(DirectMethod):
    n_crit = 10
    tauleap_threshold = 10
    micro_steps = 100

    def __init__(self, process, state, epsilon=0.03, t=0., tmax=float('inf'), steps=None, seed=None):
        super(CaoMethod, self).__init__(process, state, t=t, tmax=tmax, steps=steps, seed=seed)
        self.num_reactions = 0
        self.epsilon = epsilon
        self.rng2 = RandomState(seed)
        self.abandon_tauleap = -1

    def __iter__(self):
        while not self.has_reached_end():
            # step 1: partition reactions -- Eq. (10)
            Jcrit, Jncr = self.identify_critical_reactions()

            Irs = set(s for trans in self.transitions for s in trans.reactants)
            Incr = set(s for trans in Jncr for s in trans.reactants)

            if Incr:
                # step 2: determine noncritical tau -- Eqs. (32) and (33)
                mu = {s: sum(trans.stoichiometry.get(s, 0)*a
                             for trans, a in Jncr.items()) for s in Incr}
                var = {s: sum(trans.stoichiometry.get(s, 0)**2*a
                              for trans, a in Jncr.items()) for s in Incr}
                eps = {s: max(self.epsilon*self.state[s]*self.gi(s), 1.) for s in Incr}

                tau_ncr1 = min((eps[s]/abs(mu[s])) if mu[s] else float('inf') for s in Incr)
                tau_ncr2 = min((eps[s]**2/abs(var[s])) if var[s] else float('inf') for s in Incr)
                tau_ncr = min((tau_ncr1, tau_ncr2, self.tmax-self.time))
            else:
                tau_ncr = float('inf')

            a0 = sum(mult*prop for _, prop, mult in self.propensities.items())

            if not a0:
                break

            while True:
                # step 3: abandon tau leaping if not enough expected gain
                if tau_ncr <= self.tauleap_threshold / a0:
                    if self.abandon_tauleap == -1:
                        self.abandon_tauleap = self.step
                    it = DirectMethod.__iter__(self)
                    for _ in range(self.micro_steps):
                        trans = next(it)
                        self.num_reactions += 1
                        yield self.time, self.state, {trans: 1}
                    break
                elif self.abandon_tauleap != -1:
                    logging.debug("Abandoned tau-leaping for %d steps" % (self.step-self.abandon_tauleap))
                    self.abandon_tauleap = -1

                # step 4: determine critical tau
                ac = sum(propensity for trans, propensity in Jcrit.items())
                tau_crit = -log(self.rng.random())/ac if ac else float('inf')

                # step 5: determine actual tau
                tau = min((tau_ncr, tau_crit, self.tmax-self.time))

                # step 5a
                firings = {trans: self.rng2.poisson(propensity*tau)
                           for trans, propensity in Jncr.items()}
                firings = {trans: n for trans, n in firings.items() if n}

                # step 5b
                if tau == tau_crit:
                    # fire exactly one ciritical reaction
                    transition = None
                    pick = self.rng.random()*ac
                    for transition, propensity in Jcrit.items():
                        pick -= propensity
                        if pick < 0.:
                            break
                    firings[transition] = 1

                new_reactions = sum(firings.values())

                # avoid overshooting self.steps
                if self.steps and self.step+new_reactions > self.steps:
                    tau_ncr /= 2
                    continue

                all_reactants = sum((n*trans.true_reactants
                                     for trans, n in firings.items()), multiset())
                all_products = sum((n*trans.true_products
                                    for trans, n in firings.items()), multiset())
                net_reactants = all_reactants - all_products
                net_products = all_products - all_reactants

                # step 6a: avoid negative copy numbers
                if any(self.state[s]<n
                       for s, n in net_reactants.items()):
                    tau_ncr /= 2
                    continue

                else:
                    # step 6b: perform transitions
                    self.time += tau
                    self.step += 1
                    self.num_reactions += new_reactions
                    self.state -= net_reactants
                    for rule in self.process.rules:
                        for trans in rule.infer_transitions(net_products, self.state):
                            trans.rule = rule
                            self.add_transition(trans)
                    self.state += net_products

                    # update propensities
                    affected_species = set().union(*(trans.affected_species
                                                     for trans in firings))
                    affected = self.dependency_graph.affected_transitions(affected_species)
                    self.update_propensities(affected)
                    self.prune_transitions()

                    yield self.time, self.state, firings
                    break

        if self.step != self.steps and self.tmax < float('inf'):
            self.time = self.tmax

    def identify_critical_reactions(self):
        def critical(trans, propensity):
            if not propensity:
                return False
            elif not trans.true_reactants:
                return False
            elif L(trans) < self.n_crit:
                return True
            elif any(self.state[s] < rule.order for rule in self.process.rules for s in trans.true_products):
                return True
            else:
                return False

        def L(trans):
            return min(-float(self.state[s])/trans.stoichiometry.get(s, 0)
                       for s in trans.true_reactants
                       if trans.stoichiometry.get(s, 0) < 0)

        crit = {}
        noncrit = {}
        for trans, a, mult in self.propensities.items():
            partition = crit if critical(trans, a) else noncrit
            partition[trans] = mult*a
        return crit, noncrit

    def gi(self, species):
        g = 0
        for trans in self.transitions:
            if species not in trans.reactants:
                continue
            order = sum(trans.reactants.values())
            if order == 1:
                g = max(g, 1)
            elif order == 2 and trans.reactants[species] == 1:
                g = max(g, 2)
            elif order == 2 and trans.reactants[species] == 2:
                g = max(g, 2+(self.state[species]-1)**-1) if self.state[species] >= 2 else 3
            elif order == 3 and trans.reactants[species] == 1:
                g = max(g, 3)
            elif order == 3 and trans.reactants[species] == 2:
                g = max(g, 1.5*(2+(self.state[species]-1)**-1)) if self.state[species] >= 2 else 4.5
            elif order == 3 and trans.reactants[species] == 3:
                g = max(g, 3+(self.state[species]-1)**-1+(self.state[species]-2)**-1) if self.state[species] >= 3 else 5.5
            else:
                raise RuntimeError("Tau-leaping not implemented for reactions of order %d" % order)
        return g


if sys.argv[0].endswith('stocal/examples/validation.py'):
    # Provide some specialized methods with fixed epsilon values

    class CaoMethod_003(CaoMethod):
        def __init__(self, process, state, t=0., tmax=float('inf'), steps=None, seed=None):
            super(CaoMethod_003, self).__init__(process, state, 0.03, t=t, tmax=tmax, steps=steps, seed=seed)
    
    class CaoMethod_001(CaoMethod):
        def __init__(self, process, state, t=0., tmax=float('inf'), steps=None, seed=None):
            super(CaoMethod_001, self).__init__(process, state, 0.01, t=t, tmax=tmax, steps=steps, seed=seed)
    
    class CaoMethod_0003(CaoMethod):
        def __init__(self, process, state, t=0., tmax=float('inf'), steps=None, seed=None):
            super(CaoMethod_0003, self).__init__(process, state, 0.003, t=t, tmax=tmax, steps=steps, seed=seed)
    
    class CaoMethod_0001(CaoMethod):
        def __init__(self, process, state, t=0., tmax=float('inf'), steps=None, seed=None):
            super(CaoMethod_0001, self).__init__(process, state, 0.001, t=t, tmax=tmax, steps=steps, seed=seed)
    
    class CaoMethod_00003(CaoMethod):
        def __init__(self, process, state, t=0., tmax=float('inf'), steps=None, seed=None):
            super(CaoMethod_00003, self).__init__(process, state, 0.0003, t=t, tmax=tmax, steps=steps, seed=seed)
    
    class CaoMethod_00001(CaoMethod):
        def __init__(self, process, state, t=0., tmax=float('inf'), steps=None, seed=None):
            super(CaoMethod_00001, self).__init__(process, state, 0.0001, t=t, tmax=tmax, steps=steps, seed=seed)


# testing
import unittest
from stocal.tests.test_algorithms import TestTrajectorySampler

class TestCaoMethod(TestTrajectorySampler):
    """Test CaoMethod

    This tests the regular StochasticSimulationAlgorithm interface."""
    Sampler = CaoMethod

    @unittest.skip("Sampler does not adhere to specification")
    def test_add_transition_enables_transition(self):
        self.fail("CaoMethod violates current StochasticSimulationAlgorithm specification.")

    @unittest.skip("Sampler does not adhere to specification")
    def test_update_state_enables_infered(self):
        self.fail("CaoMethod violates current StochasticSimulationAlgorithm specification.")

    @unittest.skip("Sampler does not adhere to specification")
    def test_update_state_enables_static(self):
        self.fail("CaoMethod violates current StochasticSimulationAlgorithm specification.")


if __name__ == '__main__':
    unittest.main()
