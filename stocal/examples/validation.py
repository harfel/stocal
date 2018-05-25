"""This python script performs validation of the stocal suite

The module is meant to be run from command line as

> python stocal/examples/validation.py {run|report}

see python stocal/examples/validation.py -h

for more information.
"""
import warnings

class DataStore(object):
    def __init__(self, path):
        pass

    def add_result(self, config, result):
        warnings.warn("XXX results are not recorded yet.", RuntimeWarning)


def import_by_name(name):
    """import and return a module member given by name
    
    e.g. 'stocal.algorithms.DirectMethod' will return the class
    <class 'stocal.algorithms.DirectMethod'>
    """
    components = name.split('.')
    mod = __import__(components[0])
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod


def get_implementations(module, cls):
    return [
        member for member in module.__dict__.values()
        if isclass(member)
        and issubclass(member, cls)
        and not isabstract(member)
    ]


def run_simulation(model, algorithm):
    """Simulate model.process with algorithm

    The model.process is sampled every time unit for 50 time units.
    
    """
    def sample(trajectory, species, dt=0):
        """Sample species along trajectory every dt time units

        species is a list of species labels that should be sampled.

        Returns a tuple of two elements, the first is a list of all firing
        times, the second a dictionary that holds for each species the
        list of copy numbers at each corresponding time point. If dt is
        given, it specifies the interval at which the trajectory is sampled.
        """
        def every(trajectory, dt):
            tmax = trajectory.tmax
            trajectory.tmax = trajectory.time
            while trajectory.time < tmax:
                transitions = {}
                if trajectory.steps and trajectory.step >= trajectory.steps:
                    break
                trajectory.tmax += dt
                for trans in trajectory:
                    transitions[trans] = transitions.get(trans, 0) + 1
                yield transitions

        times = [trajectory.time]
        numbers = {s:[trajectory.state[s]] for s in species}
        it = every(trajectory, dt) if dt else iter(trajectory)
        for _ in it:
            times.append(trajectory.time)
            for s in species:
                numbers[s].append(trajectory.state[s])
        return times, numbers

    trajectory = algorithm(model.process, model.initial_state, tmax=50)
    return sample(trajectory, model.species, dt=1)


def run_validation(args):
    # initialize validation data directory/file
    args.store
    # [...]

    # collect algorithms to validate
    if not args.algo:
        from stocal import algorithms
        args.algo = get_implementations(algorithms, algorithms.TrajectorySampler)

    # collect models for validation
    if not args.models:
        from stocal.examples import dsmts
        args.models = get_implementations(dsmts, dsmts.DSMTS_Test)

    # generate N run configurations
    # from all cartesian crossings of args.models with args.algo
    from itertools import product, cycle, islice
    configurations = islice(cycle(product(args.algo, args.models)), args.N)


    # simulate in worker pool
    for algo, model in configurations:
        # XXX progress logging
        result = run_simulation(model, algo)
        # print model, algo, result
        args.store.add_result(result, (model, algo))


def report_validation(args):
    raise RuntimeError("XXX not implemented yet.")


if __name__ == '__main__':
    import sys
    import argparse
    from inspect import isclass, isabstract

    parser = argparse.ArgumentParser(prog=sys.argv[0])
    # global options
    parser.add_argument('--dir', dest='store',
                        type=DataStore,
                        default=DataStore(''),
                        help='directory for/with simulation results')
    parser.add_argument('--algo',
                        type=import_by_name,
                        action='append',
                        help='specify algorithm to be validated')
    parser.add_argument('--model',
                        action='append',
                        dest='models',
                        help='specify model to be validated')
    parser.add_argument('--cpu', metavar='N',
                        type=int,
                        default=1,
                        help='number of parallel processes')
    subparsers = parser.add_subparsers(help='validation sub-command')

    # parser for the "run" command
    parser_run = subparsers.add_parser('run', help='run simulations to generate validation data')
    parser_run.add_argument('N',
                            type=int,
                            help='number of simulations to be performed in total')
    parser_run.set_defaults(func=run_validation)

    # parser for the "report" command
    parser_report = subparsers.add_parser('report', help='assemble report from generated data')
    parser_report.set_defaults(func=report_validation)

    # parse and act
    args = parser.parse_args()
    args.func(args)
