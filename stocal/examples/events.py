import stocal


process = stocal.Process([
	stocal.MassAction(['A', 'A'], ['A2'], 0.01),
	stocal.MassAction(['A2'], ['A', 'A'], 1.),
	stocal.Event([], ['A'], 0., 1.)
])


if __name__ == '__main__' :
	traj = process.trajectory({}, tmax=100)
	for _ in traj :
		print traj.time, traj.state.get('A',0), traj.state.get('A2',0)
