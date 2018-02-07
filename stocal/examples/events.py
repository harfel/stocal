"""Event example

stocal.Event's can be added to a processes definition just like
Reactions. Process.trajectory returns an TrajectorySampler that
can cope with deterministic transitions (e.g. FirstReactionMethod).
Sampler selection and usage is entirely transparent to the user.
"""
import stocal


process = stocal.Process([
    stocal.MassAction(['A', 'A'], ['A2'], 0.01),
    stocal.MassAction(['A2'], ['A', 'A'], 1.),
    stocal.Event([], ['A'], 0., 1.)
])


if __name__ == '__main__':
    traj = process.trajectory({}, tmax=100)
    for _ in traj:
        print(traj.time, traj.state['A'], traj.state['A2'])
