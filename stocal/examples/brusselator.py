"""Brusselator

A stochastic realization of the famous Brusselator system, first
proposed in I. Prigogine and R. Lefever, Symmetry Breaking Instabilities
in Dissipative Systems, J. Chem. Phys. 48, 1695 (1968).

This is a simple example of a process with only static (non-infered)
reactions. The deterministic system exhibits ascillations when b>a+1.
"""

import stocal

a = 2.
b = 10.

process = stocal.Process([
    stocal.MassAction({}, {"x": 1}, a),
    stocal.MassAction({"x": 2, "y": 1}, {"x": 3}, 1.),
    stocal.MassAction({"x": 1}, {"y": 1, "c": 1}, b),
    stocal.MassAction({"x": 1}, {"d": 1}, 1.),
])

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    traj = process.trajectory({}, tmax=50)
    for s in traj.species:
        plt.step(traj.times, traj['y'])
    plt.show()
