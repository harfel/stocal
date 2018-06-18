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
                    try:
                        config = pickle.load(open(fname)).config
                        yield fname, self.get_stats(config)
                    except Exception as exc:
                        logging.warn("Could not access data in %s" % fname)
                        logging.info(exc, exc_info=True)
                        yield fname, None

    def get_path_for_config(self, config):
        model, algo = config
        prefix = '-'.join((algo.__name__, model.__name__))
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


def run_simulation(config):
    # setup model and algorithm
    Model, Algorithm = config
    model = Model()
    trajectory = Algorithm(model.process, model.initial_state,
                           tmax=model.tmax)

    # perform simulation
    logging.debug("Start simulation of %s with %s."
                  % (Model.__name__, Algorithm.__name__))
    result = model(trajectory)
    logging.debug("Simulation of %s with %s finished."
                  % (Model.__name__, Algorithm.__name__))
    return result


def run_in_process(queue, locks, store):
    while True :
        config = queue.get()
        if not config: break

        result = run_simulation(config)
        
        with locks[config]:
            store.feed_result(result, config)

    logging.debug("Worker finished")


def run_validation(args):
    from multiprocessing import Process, Queue, Lock
    from itertools import product

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
    def configurations(N):
        import random
        required = {
            config: (max(0, N-args.store.get_stats(config).runs)
                    if os.path.exists(args.store.get_path_for_config(config))
                    else N)
            for config in product(args.models, args.algo)
        }
        while required:
            config = random.choice(required.keys())
            required[config] -= 1
            if not required[config]:
                del required[config]
            yield config

    if args.cpu > 1:
        queue = Queue(maxsize=args.cpu)
        locks = {
            config: Lock()
            for config in product(args.models, args.algo)
        }
        processes = [Process(target=run_in_process,
                             args=(queue, locks, args.store))
                     for _ in range(args.cpu) ]
        for proc in processes:
            proc.start()
        logging.debug("%d processes started." % args.cpu)
        for config in configurations(args.N): # XXX can raise EOFError
            queue.put(config)
        logging.debug("All jobs requested.")
        for _ in processes:
            queue.put(None)
            logging.debug("Shutdown signal sent.")
        queue.close()
        for proc in processes:
            proc.join()
    else:
        for config in configurations(args.N):
            result = run_simulation(config)
            args.store.feed_result(result, config)
    logging.debug("Done.")


def report_validation(args):
    for fname, stats in args.store:
        if stats:
            figname = fname[:-len('.dat')]+'.'+args.frmt
            if not os.path.exists(figname) or os.path.getmtime(fname)>os.path.getmtime(figname):
                # only if .dat newer than 
                logging.info("Generate figure for %s" % fname)
                generate_figure(stats, figname)


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
    colormaps = [plt.cm.winter, plt.cm.copper]

    ax = fig.add_subplot(131)
    plt.title("simulation results")
    for (species, mu), cm in zip(stats.mean.items(), colormaps):
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
    plt.legend(loc=0)

    ax = fig.add_subplot(132)
    plt.title("convergence toward mean")
    ax.set_xscale('log')
    ax.set_yscale('log')
    for s, cm in zip(stats.mean, colormaps):
        ax.set_prop_cycle(plt.cycler('color', (cm(x/stats.times[-1]) for x in stats.times)))
        for t in stats.times:
            exp = rep_means[s][int(t)]
            ys = [abs((sim-exp)/rep_stdevs[s][int(t)]**2) if rep_stdevs[s][int(t)] else 0 for sim in stats.conv_mean[s, t]]
            if not all(ys): continue
            if t in stats.times[:2] or t == max(stats.times):
                ax.plot(Ns, ys, alpha=0.67, label="time = %.1f" % t)
            else:
                ax.plot(Ns, ys, alpha=0.67)
    ymin, ymax = plt.gca().get_ylim()
    ax.plot(Ns, [3/sqrt(n) for n in Ns], color='r', label="bound")
    plt.xlabel('samples N')
    plt.ylabel('normalized error')
    plt.legend(loc=3)

    ax = fig.add_subplot(133)
    plt.title("convergence toward std. dev.")
    ax.set_xscale('log')
    ax.set_yscale('log')
    for s, cm in zip(stats.mean, colormaps):
        ax.set_prop_cycle(plt.cycler('color', (cm(x/stats.times[-1]) for x in stats.times)))
        for t in stats.times:
            exp = rep_stdevs[s][int(t)]
            ys = [abs(sim/exp**2-1) if exp else 0 for sim in stats.conv_stdev[s, t]]
            if not all(ys): continue
            if t in stats.times[:2] or t == max(stats.times):
                ax.plot(Ns, ys, alpha=0.67, label="time = %.1f" % t)
            else:
                ax.plot(Ns, ys, alpha=0.67)
    ymin, ymax = plt.gca().get_ylim()
    ax.plot(Ns, [5/sqrt(n/2.) for n in Ns], color='r', label="bound")
    plt.xlabel('samples N')
    plt.ylabel('normalized error')
    plt.legend(loc=3)

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
    default_store = os.path.join('validation_data', git_label)

    parser = argparse.ArgumentParser(prog=sys.argv[0])
    # global options
    parser.add_argument('--dir', dest='store',
                        type=DataStore,
                        default=DataStore(default_store),
                        help='directory for/with simulation results')
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
    logging.basicConfig(level=logging.DEBUG)
    args = parser.parse_args()
    args.func(args)
