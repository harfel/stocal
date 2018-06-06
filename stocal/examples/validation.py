"""This python script performs validation of the stocal suite

The module is meant to be run from command line as

> python stocal/examples/validation.py {run|report}

see python stocal/examples/validation.py -h

for more information.
"""
import sys
import os
import warnings
import logging

from collections import namedtuple
from math import sqrt

try:
    import numpy as np
except ImportError:
        logging.error("Example validation.py requires numpy.")
        sys.exit(1)


class Stats(namedtuple(
    '_Stats',
    ('runs', 'times', 'mean', 'M2', 'conv_mean', 'conv_stdev', 'config'))):

    @property
    def stdev(self):
        return {
            s: (values/(self.runs-1))**.5
            for s, values in self.M2.items()
        }


class DataStore(object):
    checkpoints = [int(sqrt(10)**n) for n in range(100)][1:]

    def __init__(self, path):
        import errno
        self.path = path
        try:
            os.makedirs(self.path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    def __iter__(self):
        import pickle
        for dirpath, _, filenames in os.walk(self.path):
            for name in filenames:
                fname = os.path.join(dirpath, name)
                if fname.endswith('.dat'):
                    config = pickle.load(open(fname)).config
                    yield fname, self.get_stats(config)

    def get_path_for_config(self, config):
        model, algo = config
        prefix = '-'.join((model.__name__, algo.__name__))
        return os.path.join(self.path, prefix+'.dat')

    def feed_result(self, result, config):
        # online aggregation adapted from
        # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm

        import pickle
        from shutil import copyfile
        from math import log, floor

        fname = self.get_path_for_config(config)
        if os.path.exists(fname):
            stats = pickle.load(open(fname))
            N = stats.runs + 1
            times = stats.times
            delta = {
                s: values - stats.mean[s]
                for s, values in result[1].items()
            }
            mean = {
                s: values + delta[s]/float(N)
                for s, values in stats.mean.items()
            }
            delta2 = {
                s: values - mean[s]
                for s, values in result[1].items()
            }
            M2 = {
                s: values + delta[s]*delta2[s]
                for s, values in stats.M2.items()
            }
            conv_mean = stats.conv_mean
            conv_stdev = stats.conv_stdev
            copyfile(fname, fname+'~')

        else:
            N = 1
            times, mean = result
            M2 = {s: np.array([0. for _ in mean[s]]) for s in mean}
            conv_mean = {}
            conv_stdev = {}

        if N in self.checkpoints:
            for s, means in mean.items():
                for t, mu, m2 in zip(times, means, M2[s]):
                    if (s, t) not in conv_mean:
                        conv_mean[s, t] = np.array([])
                        conv_stdev[s, t] = np.array([])
                    conv_mean[s, t] = np.append(conv_mean[s, t], [mu])
                    conv_stdev[s, t] = np.append(conv_stdev[s, t], [sqrt(m2/(N-1))])

        with open(fname, 'w') as outfile:
            stats = Stats(N, times, mean, M2, conv_mean, conv_stdev, config)
            outfile.write(pickle.dumps(stats))

        if os.path.exists(fname+'~'):
            os.remove(fname+'~')

    def get_stats(self, config):
        import pickle
        model, algo = config
        fname = self.get_path_for_config(config)
        return pickle.load(open(fname))


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
        numbers = {s: np.array([trajectory.state[s]]) for s in species}
        it = every(trajectory, dt) if dt else iter(trajectory)
        for _ in it:
            times.append(trajectory.time)
            for s in species:
                numbers[s] = np.append(numbers[s], [trajectory.state[s]])
        return np.array(times), numbers

    trajectory = algorithm(model.process, model.initial_state, tmax=50)
    return sample(trajectory, model.species, dt=1)


def run_validation(args):
    def get_implementations(module, cls):
        return [
            member for member in module.__dict__.values()
            if isclass(member)
            and issubclass(member, cls)
            and not isabstract(member)
        ]

    # collect algorithms to validate
    if not args.algo:
        from stocal import algorithms
        args.algo = get_implementations(algorithms, algorithms.TrajectorySampler)

    # collect models for validation
    if not args.models:
        from stocal.examples.dsmts import models as dsmts
        args.models = get_implementations(dsmts, dsmts.DSMTS_Test)

    # generate required simulation configurations
    def required_runs(config, N):
        try:
            return max(0, args.store.get_stats(config).runs-N)
        except IOError:
            return N
    from itertools import product, repeat, chain
    configurations = chain(*(
        repeat(config, required_runs(config, args.N))
        for config in product(args.models, args.algo)
    ))

    # simulate in worker pool
    for model, algo in configurations:
        result = run_simulation(model, algo)
        args.store.feed_result(result, (model, algo))


def report_validation(args):
    for fname, stats in args.store:
        generate_figure(stats, fname[:-len('.dat')]+'.'+args.frmt)


def generate_figure(stats, fname):
    try:
        from matplotlib import pyplot as plt
    except ImportError:
        logging.error("Example validation.py requires matplotlib.")
        sys.exit(1)

    model, algo = stats.config
    rep_times, rep_means = model.reported_means()
    rep_times, rep_stdevs = model.reported_stdevs()
    Ns = DataStore.checkpoints[:len(stats.conv_mean.values()[0])]

    fig = plt.figure(figsize=plt.figaspect(.3))
    title = '%s %s (%d samples)' % (model.__name__, algo.__name__, stats.runs)
    fig.suptitle(title)
    cm = plt.cm.winter

    ax = fig.add_subplot(131)
    plt.title("simulation results")
    for species, mu in stats.mean.items():
        low = mu - stats.stdev[species]**.5
        high = mu + stats.stdev[species]**.5

        rep_low = rep_means[species] - rep_stdevs[species]
        rep_high = rep_means[species] + rep_stdevs[species]

        ax.fill_between(stats.times, rep_low, rep_high, facecolor=cm(0), alpha=0.3)
        ax.fill_between(stats.times, low, high, facecolor=cm(0.99), alpha=0.3)
        ax.plot(rep_times, rep_means[species], color=cm(0), label='$\mathregular{%s_{exp}}$' % species)
        ax.plot(stats.times, mu, color=cm(0.99), label='$\mathregular{%s_{sim}}$' % species, alpha=.67)
    plt.xlabel('time')
    plt.ylabel('# molecules')
    plt.legend()

    ax = fig.add_subplot(132)
    plt.title("convergence toward mean")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_prop_cycle(plt.cycler('color', (cm(x/stats.times[-1]) for x in stats.times)))
    for s in stats.mean:
        for t in stats.times:
            exp = rep_means[s][int(t)]
            ys = [abs((sim-exp)/rep_stdevs[s][int(t)]**2) if rep_stdevs[s][int(t)] else 0 for sim in stats.conv_mean[s, t]]
            if not all(ys): continue
            ax.plot(Ns, ys, alpha=0.67)
    ymin, ymax = plt.gca().get_ylim()
    ax.plot(Ns, [3/sqrt(n) for n in Ns], color='r')
    plt.xlabel('samples N')
    plt.ylabel('normalized error')
    plt.legend()

    ax = fig.add_subplot(133)
    plt.title("convergence toward std. dev.")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_prop_cycle(plt.cycler('color', (cm(x/stats.times[-1]) for x in stats.times)))
    for s in stats.mean:
        for t in stats.times:
            exp = rep_stdevs[s][int(t)]
            ys = [abs(sim/exp**2-1) if exp else 0 for sim in stats.conv_stdev[s, t]]
            if not all(ys): continue
            ax.plot(Ns, ys, alpha=0.67)
    ymin, ymax = plt.gca().get_ylim()
    ax.plot(Ns, [5/sqrt(n/2.) for n in Ns], color='r')
    plt.xlabel('samples N')
    plt.ylabel('normalized error')
    plt.legend()

    fig.savefig(fname)


if __name__ == '__main__':
    import sys
    import argparse
    import subprocess
    from inspect import isclass, isabstract

    def import_by_name(name):
        """import and return a module member given by name
        
        e.g. 'stocal.algorithms.DirectMethod' will return the class
        <class 'stocal.algorithms.DirectMethod'>
        """
        from importlib import import_module
        module, member = name.rsplit('.', 1)
        mod = import_module(module)
        return getattr(mod, member)

    git_label = subprocess.check_output(["git", "describe"]).strip()
    default_store = os.path.join('validation', git_label)

    parser = argparse.ArgumentParser(prog=sys.argv[0])
    # global options
    parser.add_argument('--dir', dest='store',
                        type=DataStore,
                        default=DataStore(default_store),
                        help='directory for/with simulation results')
    #parser.add_argument('--cpu', metavar='N',
    #                    type=int,
    #                    default=1,
    #                    help='number of parallel processes')
    
    subparsers = parser.add_subparsers(help='validation sub-command')

    # parser for the "run" command
    parser_run = subparsers.add_parser('run', help='run simulations to generate validation data')
    parser_run.add_argument('N',
                            type=int,
                            help='number of simulations to be performed in total')
    parser_run.add_argument('--algo',
                            type=import_by_name,
                            action='append',
                            help='specify algorithm to be validated')
    parser_run.add_argument('--model',
                            type=import_by_name,
                            action='append',
                            dest='models',
                            help='specify model to be validated')
    parser_run.set_defaults(func=run_validation)

    # parser for the "report" command
    parser_report = subparsers.add_parser('report', help='generate figures from generated data')
    parser_report.add_argument('--format',
                               action='store',
                               dest='frmt',
                               default='png',
                               help='file format of generated figures')
    parser_report.set_defaults(func=report_validation)

    # parse and act
    args = parser.parse_args()
    args.func(args)
