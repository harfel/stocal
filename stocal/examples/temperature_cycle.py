"""Example of a process with time-dependent reactions rates

This example simulates the association and dissociation of two molecules
X into a dimer X2. We assume regular mass action kinetics for the
association, but model the dissociation to be temperature dependent.
To do so, we observe that the quilibrium constant

    K = k_forward/k_backward = exp(-(dH-T*dS)/T)

is temperature dependent and allows us calculate a temperature-dependent
backward rate for a fixed forward rate.
"""
from math import exp, log, sin, pi
import stocal

k_forward = .1      # the fixed forward rate
x_tot = 100         # total number of molecules
dG = log(x_tot/3.)  # a free energy that leads to equipartition for T=1
dH = -10            # enthalpy change of association
dS = dH - dG        # entropy change of association

def temp(time, low=0.5, high=1.5, period=50):
    """Temperature cycle"""
    return low+(high-low)*(sin(2*pi*time/period)+1)/2.


class Dissociation(stocal.MassAction):
    """Temperature dependent dissocition

    The normal mass action propensity is modulated by a factor that
    accounts for the system temperature at any time. The stochastic
    rate constant is therefore k = k_forward * exp((dH-T*dS)/T).
    k_forward is provided as rate constant to the MassAction constructor
    and the exponential transform is performed in the overloaded
    propensity method.
    """
    def propensity(self, state, time):
        T = temp(time)
        return super(Dissociation, self).propensity(state) * exp((dH-T*dS)/T)


process = stocal.Process([
    stocal.MassAction(['x', 'x'], ['x2'], 2*k_forward),
    Dissociation(['x2'], ['x', 'x'], k_forward),
])

state = {'x2': x_tot//3, 'x': x_tot-2*x_tot//3}

if __name__ == '__main__':
    traj = process.sample(state, tmax=125.)
    for trans in traj:
        print(traj.time, traj.state['x'], traj.state['x2'], temp(traj.time))
