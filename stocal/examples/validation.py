"""This python script performs validation of the stocal suite

The module is meant to be run from command line as

> python stocal/examples/validation.py {run|report}

see python stocal/examples/validation.py -h

for more information.
"""
import os
import warnings
import logging

from math import sqrt


class DataStore(object):
    checkpoints = [int(sqrt(10)**n) for n in range(100)][1:]

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
        from shutil import copyfile
        from math import log, floor

        fname = self.get_path_for_config(config)
        if os.path.exists(fname):
            N, times, mean, M2, conv_mean, conv_stdev = pickle.load(open(fname))
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
            copyfile(fname, fname+'~')

        else:
            N = 1
            times = result[0]
            mean = {
                s: [float(x) for x in values]
                for s, values in result[1].items()
            }
            M2 = {s: [0. for _ in mean[s]] for s in mean}
            conv_mean = {}
            conv_stdev = {}

        if N in self.checkpoints:
            for s, means in mean.items():
                for t, mu, m2 in zip(times, means, M2[s]):
                    if (s, t) not in conv_mean:
                        conv_mean[s, t] = []
                        conv_stdev[s, t] = []
                    conv_mean[s, t].append(mu)
                    conv_stdev[s, t].append(sqrt(m2/(N-1)))

        with open(fname, 'w') as outfile:
            outfile.write(pickle.dumps((N, times, mean, M2, conv_mean, conv_stdev)))

    def get_stats(self, config):
        import pickle
        model, algo = config
        fname = self.get_path_for_config(config)
        N, times, mean, M2, conv_mean, conv_stdev = pickle.load(open(fname))
        stdev = {
            s: [sqrt(m2/(N-1)) for m2 in values]
            for s, values in M2.items()
        }
        return N, times, mean, stdev, conv_mean, conv_stdev


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
            generate_figure(data, config, fname)


def generate_figure(data, config, fname):
    try:
        from matplotlib import pyplot as plt
    except ImportError:
        logging.error("Example validation.py requires matplotlib.")

    cm = plt.cm.winter

    model, algo = config
    N, times, mean, stdev, conv_mean, conv_stdev = data

    rep_times, rep_means = model.reported_means()
    rep_times, rep_stdevs = model.reported_stdevs()

    Ns = DataStore.checkpoints[:len(conv_mean.values()[0])]

    fig = plt.figure(figsize=plt.figaspect(.3))
    title = '%s %s (%d samples)' % (model.__name__, algo.__name__, N)
    fig.suptitle(title)

    ax = fig.add_subplot(131)
    plt.title("simulation results")
    for species, mu in mean.items():
        low = [y-sqrt(s) for y,s in zip(mu, stdev[species])]
        high = [y+sqrt(s) for y,s in zip(mu, stdev[species])]

        rep_low = [m-y for m,y in zip(rep_means[species], rep_stdevs[species])]
        rep_high = [m+y for m,y in zip(rep_means[species], rep_stdevs[species])]

        ax.fill_between(times, rep_low, rep_high, facecolor=cm(0), alpha=0.3)
        ax.fill_between(times, low, high, facecolor=cm(0.99), alpha=0.3)
        ax.plot(rep_times, rep_means[species], color=cm(0), label='$\mathregular{%s_{exp}}$' % species)
        ax.plot(times, mu, color=cm(0.99), label='$\mathregular{%s_{sim}}$' % species, alpha=.67)
    plt.xlabel('time')
    plt.ylabel('# molecules')
    plt.legend()

    ax = fig.add_subplot(132)
    plt.title("convergence toward mean")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_prop_cycle(plt.cycler('color', (cm(x/times[-1]) for x in times)))
    for s in mean:
        for t in times:
            exp = rep_means[s][int(t)]
            ys = [abs((sim-exp)/rep_stdevs[s][int(t)]**2) if rep_stdevs[s][int(t)] else 0 for N, sim in zip(Ns, conv_mean[s, t])]
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
    ax.set_prop_cycle(plt.cycler('color', (cm(x/times[-1]) for x in times)))
    for s in mean:
        for t in times:
            exp = rep_stdevs[s][int(t)]
            ys = [abs(sim/exp**2-1) if exp else 0 for N, sim in zip(Ns, conv_stdev[s, t])]
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
        from stocal.examples.dsmts import models as dsmts
        args.models = get_implementations(dsmts, dsmts.DSMTS_Test)


    args.func(args) # XXX run in worker pool
