import numpy as np
import stocal as sc
from matplotlib import pyplot as plt

process = sc.Process([
	sc.MassAction('2*A -> A2', 10.), # support string syntax to specify reactions
	sc.MassAction('A2 -> 2*A', 1.),
	sc.MassAction('-> A', 1.),
])

# Generate and plot individual trajectory
traj = process.trajectory({}, tmax=1000)
plt.plot(traj.times, traj['A'])
plt.plot(traj.times, traj['A2'])
plt.show()

# Generate homogenized trajectory by passing a sequence of times
times = np.linspace(0,1000,1001)
traj = process.trajectory({}, times)
plt.plot(times, traj['A'])
plt.show()

# Generate an ensemble
# Ensembles store information about average and standard deviation
# of a process. Ensembles retain a given number of specimen trajectories.
# The following draws the process average (bold line), ten specimen
# trajectories (faint lines), and a shaded area demarking standard
# deviations.
ensemble = process.ensemble({}, times, samples=100, specimen=10)
plt.fill_between(times, ensemble.lower['A'], ensemble.upper['A'], alpha=.67)
plt.plot(times, ensemble.specimen['A'])
plt.plot(times, ensemble.average['A'], linewidth=2)
plt.show()

# Long running calculations can be performed in the background through
# the stocal configuration. 
sc.config(background=True)


# The call to Process.ensemble will now return immediately. The returned
# proxy object will evaluate ensemble properties of all the calculations
# that have been executed up to that point.
ensemble = process.ensemble(process, {}, times, samples=1000)
plt.plot(times, ensemble.average['A'])

ensemble.task # Access to the running task (ETA, pause/resume, etc.)

# For long running simulations, persistence would be real nice:
# This can be done by declaring a stocal session. A session is simply
# a path to a directory, where persistent data can be stored.
# Sessions retain calculated trajectories, ensembles, etc. so that they
# do not need to be recalculated.
sc.config(session='./store', seed=42)
