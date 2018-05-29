"""This python script performs validation of the stocal suite

The module is meant to be run from command line as

> python stocal/examples/validation.py {run|report}

see python stocal/examples/validation.py -h

for more information.
"""
import os
import warnings
import logging

try:
    from matplotlib import pyplot
except ImportError:
    logging.error("Example validation.py requires matplotlib.")


class DataStore(object):
    def __init__(self, path):
        import subprocess
        import errno
        git_label = subprocess.check_output(["git", "describe"]).strip()
        self.dir = os.path.join(path, git_label)

        try:
            os.makedirs(self.dir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    def get_path_for_config(self, config):
        model, algo = config
        prefix = '-'.join((model.__name__, algo.__name__))
        return os.path.join(self.dir, prefix+'.dat')

    def feed_result(self, result, config):
        # online aggregation adapted from
        # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm

        import pickle
        fname = self.get_path_for_config(config)
        if os.path.exists(fname):
            N, times, mean, M2 = pickle.load(open(fname))
            assert result[0] == times
            N += 1
            # delta = result - mean
            delta = {
                s: [v-x for v, x in zip(values, mean[s])]
                for s, values in result[1].items()
            }
            # mean += delta / N
            mean = {
                s: [m+d/N for m, d in zip(mean[s], delta[s])]
                for s, values in mean.items()
            }
            # delta2 = result  - mean
            delta2 = {
                s: [x-m for x, m in zip(values, mean[s])]
                for s, values in result[1].items()
            }
            # M2 += delta * delta2
            M2 = {
                s: [m2 + d1*d2 for m2, d1, d2 in zip(values, delta[s], delta2[s])]
                for s, values in M2.items()
            }

        else:
            N = 1
            times = result[0]
            mean = {
                s: [float(x) for x in values]
                for s, values in result[1].items()
            }
            M2 = {s: [0. for _ in mean[s]] for s in mean}

        with open(fname, 'w') as outfile:
            outfile.write(pickle.dumps((N, times, mean, M2))) # XXX backup copy

    def get_stats(self, config):
        import pickle
        from math import sqrt
        model, algo = config
        fname = self.get_path_for_config(config)
        N, times, mean, M2 = pickle.load(open(fname))
        stdev = {
            s: [sqrt(m2/(N-1)) for m2 in values]
            for s, values in M2.items()
        }
        return N, times, mean, stdev


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
    # generate N run configurations
    # from all cartesian crossings of args.models with args.algo
    from itertools import product, cycle, islice
    configurations = islice(cycle(product(args.models, args.algo)), args.N)


    # simulate in worker pool
    for model, algo in configurations:
        # XXX progress logging
        result = run_simulation(model, algo)
        args.store.feed_result(result, (model, algo))


def report_validation(args):
    for algo in args.algo:
        for model in args.models:
            config = model, algo
            data = args.store.get_stats(config)
            fname = args.store.get_path_for_config(config)[:-len('.dat')] + '.png' # or '.pdf'
            try:
                generate_figure(data, config, fname)
            except ValueError:
                logging.error("Could not generate report for %s" % str(config))


def generate_figure(data, config, fname):
    from math import sqrt

    model, algo = config
    N, times, avgs, var = data

    fig = pyplot.figure()
    ax = fig.add_subplot(111)

    for species, avg in avgs.items():
        low = [y-sqrt(s) for y,s in zip(avg, var[species])]
        high = [y+sqrt(s) for y,s in zip(avg, var[species])]
        title = '%s %s (%d samples)' % (model.__name__, algo.__name__, N)

        ax.fill_between(times, low, high, alpha=0.3)
        ax.plot(times, avgs[species], label=species)
        pyplot.xlabel('time')
        pyplot.ylabel('# molecules')
        pyplot.title(title)
        pyplot.legend()
    fig.savefig(fname)


if __name__ == '__main__':
    import sys
    import argparse
    from inspect import isclass, isabstract


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


    parser = argparse.ArgumentParser(prog=sys.argv[0])
    # global options
    parser.add_argument('--dir', dest='store',
                        type=DataStore,
                        default=DataStore('validation'),
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

    # collect algorithms to validate
    if not args.algo:
        from stocal import algorithms
        args.algo = get_implementations(algorithms, algorithms.TrajectorySampler)

    # collect models for validation
    if not args.models:
        from stocal.examples import dsmts
        args.models = get_implementations(dsmts, dsmts.DSMTS_Test)


    args.func(args) # XXX run in worker pool
