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
    from numpy.random import poisson
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
        self.epsilon = epsilon

    def __iter__(self):
        while not self.has_reached_end():
            # step 1: partition reactions -- Eq. (10)
            Jcrit, Jncr = self.identify_critical_reactions()

            Irs = [s for s in self.state
                   if any(s in trans.reactants for trans in self.transitions)]

            Incr = [s for s in self.state
                    if any(s in trans.reactants for trans in Jncr)]

            if Incr:
                # step 2: determine noncritical tau -- Eqs. (32) and (33)
                mu = {s: sum(self._nu(trans, s)*a
                             for trans, a in Jncr.items()) for s in Irs}
                var = {s: sum(self._nu(trans, s)**2*a
                              for trans, a in Jncr.items()) for s in Irs}

                eps = {s: max(self.epsilon*self.state[s]*self.gi(s), 1.) for s in Incr}

                tau_ncr1 = min(eps[s]/abs(mu[s]) for s in Incr)
                tau_ncr2 = min(eps[s]**2/abs(var[s]) for s in Incr)
                tau_ncr = min((tau_ncr1, tau_ncr2))
            else:
                tau_ncr = float('inf')

            # step 3: check whether to abandon tau leaping
            a0 = sum(mult*prop for _, prop, mult in self.propensities.items())

            if not a0:
                break
            
            if tau_ncr < self.tauleap_threshold / a0:
                it = DirectMethod.__iter__(self)
                for _ in xrange(self.micro_steps):
                    yield next(it)
                    # XXX yield triple rather than transition:
                    # trans = next(it)
                    # yield self.time, self.state, {trans: 1}
                continue
            else:
                # step 4: determine critical tau
                ac = sum(propensity for trans, propensity in Jcrit.items()) if tau_ncr < float('inf') else a0
                tau_crit = -log(self.rng.random())/ac if ac else float('inf')

                # step 5: determine actual tau
                tau = min((tau_ncr, tau_crit))

                # step 5a
                firings = {trans: poisson(propensity*tau) # XXX does not use self.seed
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

            # step 6: perform transitions
            self.time += tau
            all_reactants = sum((n*trans.true_reactants
                                 for trans, n in firings.items()), multiset())
            all_products = sum((n*trans.true_products
                                for trans, n in firings.items()), multiset())

            if any(self.state[s]<n
                   for s, n in (all_reactants-all_products).items()):
                # XXX account for potentially negative species
                raise RuntimeError("negative species count encountered during simulation.")

            self.step += sum(firings.values()) # XXX seperate counts for step and num_transitions
            self.state -= all_reactants
            for rule in self.process.rules:
                for trans in rule.infer_transitions(all_products, self.state):
                    trans.rule = rule
                    self.add_transition(trans)
            self.state += all_products

            # update propensities
            affected_species = set().union(*(trans.affected_species
                                             for trans in firings))
            affected = self.dependency_graph.affected_transitions(affected_species)
            self.update_propensities(affected)

            yield firings.keys()[0]
            # XXX yield triple rather than transition:
            #yield self.time, self.state, firings

        if self.step != self.steps and self.tmax < float('inf'):
            self.time = self.tmax


    def identify_critical_reactions(self):
        def critical(trans, propensity):
            if not propensity:
                return True
            elif not trans.true_reactants:
                return False
            elif L(trans) < self.n_crit:
                return True
            elif any(self.state[s] < rule.order for rule in self.process.rules for s in trans.true_products):
                return True
            else:
                return False

        def L(trans):
            return min(-float(self.state[r])/self._nu(trans, r)
                       for r in trans.true_reactants
                       if self._nu(trans, r) < 0)

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
            # XXX generalize the below
            if order == 1:
                g = max(g, 1)
            elif order == 2 and trans.reactants[species] == 1:
                g = max(g, 2)
            elif order == 2 and trans.reactants[species] == 2:
                g = max(g, 2+(self.state[species]-1)**-1)

            elif order == 3 and trans.reactants[species] == 1:
                g = max(g, 3)
            elif order == 3 and trans.reactants[species] == 2:
                g = max(g, 1.5*(2+(self.state[species]-1)**-1))
            elif order == 3 and trans.reactants[species] == 3:
                g = max(g, 3+(self.state[species]-1)**-1+(self.state[species]-2)**-1)
            else:
                raise RuntimeError("Tau-leaping not implemented for reactions of order %d" % order)
        return g

    @staticmethod
    def _nu(trans, r):
        return trans.true_products[r] - trans.true_reactants[r]


# testing
import unittest
from stocal.tests.test_algorithms import TestTrajectorySampler

class TestCaoMethod(TestTrajectorySampler):
    """Test CaoMethod

    This tests the regular TrajectorySampler interface."""
    Sampler = CaoMethod


if __name__ == '__main__':
    #from stocal.examples.pre2017 import process, initial_state

    #traj = CaoMethod(process, initial_state, 0.03, steps=sum(initial_state.values()))
    #for time, state, transitions in traj:
    #    print time, traj.step, sum(transitions.values())
    
    unittest.main()
