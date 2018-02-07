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

traj = stocal.Process([
    stocal.MassAction({}, {"x": 1}, a),
    stocal.MassAction({"x": 2, "y": 1}, {"x": 3}, 1.),
    stocal.MassAction({"x": 1}, {"y": 1, "c": 1}, b),
    stocal.MassAction({"x": 1}, {"d": 1}, 1.),
]).trajectory({}, tmax=50)

print("# time\tx\ty\tc\td")
for transition in traj:
    print(traj.time, '\t'.join(str(traj.state[s]) for s in "xycd"))
