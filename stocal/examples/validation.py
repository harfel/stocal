"""This python script performs validation of the stocal suite

The module is meant to be run from command line as

> python stocal/examples/validation.py {run|report}

see python stocal/examples/validation.py -h

for more information.

The script can be used to run validation tests from the DSMTS suite
(by default stocal.examples.dsmts.models) on stochastic simulation
algorithm implementations (by default all algorithms in
stocal.algorithms). Samples trajectories are fed into a DataStore
that performs aggregation to estimate mean and standard deviations
for all points along the trajectory. By default, the path of the
data store is determined from the current git version (via git describe)
but can be changed via command line option. Optional multiprocessing
allows to perform simulations in parallel.

Simulation results can be visualized (as png or pdf images) using the
report command. The report command also generates a validation.tex file
that can be compiled into a report.
"""
import sys
import os
import logging

from collections import namedtuple
from math import sqrt

try:
    import numpy as np
    import jinja2
except ImportError:
    logging.error("Example validation.py requires numpy and jinja2.")
    sys.exit(1)


class Stats(namedtuple('_Stats',
                       ('runs', 'times', 'mean', 'M2', 'conv_mean', 'conv_stdev', 'config'))):
    """Simulation result statistics

    Stats collects simulation statistics for use in the DataStore,
    where each algorithm/model configuration has one associated
    Stats instance.
    
    Stats groups the number of runs, a sequence of trajectory time
    points, together with mean and M2 values for each species in the
    model. (M2 is a temporary value used to keep track of standard
    deviations. See DataStore docmentation for details) mean and M2 are
    dictionaries that map each species identifier to a sequence of
    floating point numbers.

    Stats.stdev returns a dictionary of standard deviation sequences
    for each species.
    """
    @property
    def stdev(self):
        """Dictionary of standard deviation sequences for each species"""
        return {
            s: (values/(self.runs-1))**.5
            for s, values in self.M2.items()
        }


class DataStore(object):
    """Persistent store for aggregated data
    
    This class provides a data store that maintains statistics of
    simulation results throughout multiple incarnations of the
    validation script.

    A DataStore accepts individual simulation results for given
    configurations, and allows retrieval of the aggregated statistics
    for a given configuration.
    """
    # times at which current aggregate data is memorized
    checkpoints = [int(sqrt(10)**n) for n in range(100)][1:]

    def __init__(self, path):
        import errno
        self.path = path
        try:
            os.makedirs(self.path)
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise

    def __iter__(self):
        """Access all data files and statistics in the data store"""
        import pickle
        for dirpath, _, filenames in os.walk(self.path):
            for name in filenames:
                fname = os.path.join(dirpath, name)
                if fname.endswith('.dat'):
                    try:
                        with open(fname, 'rb') as fstats:
                            config = pickle.load(fstats).config
                        yield fname, self.get_stats(config)
                    except Exception as exc:
                        logging.warn("Could not access data in %s", fname)
                        logging.info(exc, exc_info=True)
                        yield fname, None

    def get_path_for_config(self, config):
        """Retrieve path of datafile for a given configuration"""
        model, algo = config
        prefix = '-'.join((algo.__name__, model.__name__))
        return os.path.join(self.path, prefix+'.dat')

    def feed_result(self, result, config):
        """Add a single simulation result for a given configuration

        feed_result uses an online algorithm to update mean and
        standard deviation with every new result fed into the store.
        (The online aggregation is adapted from
        https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm)
        
        At times defined in self.checkpoints, a dump of the current
        statistics is memorized in Stats.conv_mean and Stats.conv_stdev.
        """

        import pickle
        from shutil import copyfile
        from math import log, floor

        fname = self.get_path_for_config(config)
        if os.path.exists(fname):
            with open(fname, 'rb') as fstats:
                stats = pickle.load(fstats)
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

        with open(fname, 'wb') as outfile:
            stats = Stats(N, times, mean, M2, conv_mean, conv_stdev, config)
            outfile.write(pickle.dumps(stats))

        if os.path.exists(fname+'~'):
            os.remove(fname+'~')

    def get_stats(self, config):
        """Read stats for a given configuration"""
        import pickle
        fname = self.get_path_for_config(config)
        with open(fname, 'rb') as fstats:
            return pickle.load(fstats)


def run_simulation(Model, Algorithm, max_steps=100000):
    """Perform single simulation of Model using Algorithm.

    Returns the result of a single simulation run.
    """
    # setup model and algorithm
    model = Model()
    trajectory = Algorithm(model.process, model.initial_state,
                           tmax=model.tmax, steps=max_steps)

    # perform simulation
    logging.debug("Start simulation of %s with %s.",
                  Model.__name__, Algorithm.__name__)
    result = model(trajectory)
    logging.debug("Simulation of %s with %s finished.",
                  Model.__name__, Algorithm.__name__)
    return result


def run_in_process(queue, locks, store):
    """Worker process for parallel execution of simulations.
    
    The worker continuously fetches a simulation configuration from
    the queue, runs the simulation and feeds the simulation result
    into the data store. The worker stops if it fetches a single None
    from the queue.
    """
    while True:
        config = queue.get()
        if not config:
            break

        try:
            result = run_simulation(*config)
        except Exception as exc:
            logging.warning("Could not run simulation for %s", str(config))
            logging.info(exc, exc_info=True)

        with locks[config]:
            try:
                store.feed_result(result, config)
            except Exception as exc:
                logging.warning("Could not store result for %s", str(config))
                logging.info(exc, exc_info=True)

    logging.debug("Worker finished")


def run_validation(args):
    """Perform validation simulations.
    
    Run simulations required for the store to hold aggregregated
    statistics from args.N samples for each given algorithm and model
    combination. If args.model is not given, models classes are
    loaded from stocal.examples.dsmts.models. If args.algo is not given,
    algorithms are loaded from stocal.algorithms.

    If args.cpu is given and greater than 1, simulations are performed
    in parallel.
    """
    from multiprocessing import Process, Queue, Lock
    from inspect import isclass, isabstract
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
        args.algo = get_implementations(algorithms, algorithms.StochasticSimulationAlgorithm)

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
        configs = list(required)
        while configs:
            next_config = random.choice(configs)
            required[next_config] -= 1
            if not required[next_config]:
                configs.remove(next_config)
            yield next_config

    if args.cpu > 1:
        queue = Queue(maxsize=args.cpu)
        locks = {
            config: Lock()
            for config in product(args.models, args.algo)
        }
        processes = [Process(target=run_in_process,
                             args=(queue, locks, args.store))
                     for _ in range(args.cpu)]
        for proc in processes:
            proc.start()
        logging.debug("%d processes started." % args.cpu)
        for config in configurations(args.N):
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
            result = run_simulation(*config)
            args.store.feed_result(result, config)
    logging.debug("Done.")


def report_validation(args, frmt='png', template='doc/validation.tex'):
    """Generate figures for all results in args.store"""
    def camel_case_split(identifier):
        import re
        matches = re.finditer('.+?(?:(?<=[a-z])(?=[A-Z])|(?<=[A-Z])(?=[A-Z][a-z])|$)', identifier)
        return ' '.join(m.group(0) for m in matches)

    figures = {}

    # generate figures for the entire data store
    for fname, stats in args.store:
        if stats:
            figname = fname[:-len('.dat')]+'.'+ frmt
            if not os.path.exists(figname) or os.path.getmtime(fname) > os.path.getmtime(figname):
                # only if .dat newer than
                logging.debug("Generate figure for %s", fname)
                try:
                    generate_figure(stats, figname)
                except Exception as exc:
                    logging.warning("Could not generate figure for %s", fname)
                    logging.info(exc, exc_info=True)

            model, algo = stats.config
            algo_name = camel_case_split(algo.__name__).replace('_', ' ')
            figures[algo_name] = sorted(figures.get(algo_name, [])+[figname])

    # populate latex template
    latex_jinja_env = jinja2.Environment(
        block_start_string = '\BLOCK{',
        block_end_string = '}',
        variable_start_string = '\VAR{',
        variable_end_string = '}',
        comment_start_string = '\#{',
        comment_end_string = '}',
        line_statement_prefix = '%%',
        line_comment_prefix = '%#',
        trim_blocks = True,
        autoescape = False,
        loader = jinja2.FileSystemLoader(os.path.abspath('.'))
    )
    template = latex_jinja_env.get_template(template)
    context = {
        'version': os.path.basename(args.store.path).replace('_', '\_'),
        'methods': figures,
    }
    reportfile = args.reportfile or 'validation-%s.tex' % context['version']
    with open(reportfile, 'w') as report:
        report.write(template.render(**context))

def generate_figure(stats, fname):
    """Generate figure for given stats and save it to fname."""
    try:
        from matplotlib import pyplot as plt
    except ImportError:
        logging.error("Example validation.py requires matplotlib.")
        sys.exit(1)

    model, algo = stats.config
    rep_times, rep_means = model.reported_means()
    rep_times, rep_stdevs = model.reported_stdevs()
    Ns = DataStore.checkpoints[:len(list(stats.conv_mean.values())[0])]

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
        ax.plot(rep_times, rep_means[species],
                color=cm(0), label=r'$\mathregular{%s_{exp}}$' % species)
        ax.plot(stats.times, mu,
                color=cm(0.99), label=r'$\mathregular{%s_{sim}}$' % species, alpha=.67)
    plt.xlabel('time')
    plt.ylabel('# molecules')
    plt.legend(loc=0)

    ax = fig.add_subplot(132)
    plt.title("convergence toward mean")
    ax.set_xscale('log')
    ax.set_yscale('log')
    for s, cm in zip(stats.mean, colormaps):
        ax.set_prop_cycle(plt.cycler('color',
                                     (cm(x/stats.times[-1]) for x in stats.times)))
        for t in stats.times:
            exp = rep_means[s][int(t)]
            ys = [abs((sim-exp)/rep_stdevs[s][int(t)]**2)
                  if rep_stdevs[s][int(t)]
                  else 0
                  for sim in stats.conv_mean[s, t]]
            if not all(ys):
                continue
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
            if not all(ys):
                continue
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
    plt.close()


if __name__ == '__main__':
    import argparse
    import subprocess

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

    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        description="Perform statistical validation tests.",
        epilog="""If --dir is provided, it specifies a directory used to
               hold validation data.""")
    # global options
    parser.add_argument('--dir', dest='store',
                        type=DataStore,
                        default=DataStore(default_store),
                        help='directory for/with simulation results')

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
    parser_run.add_argument('--cpu', metavar='N',
                            type=int,
                            default=1,
                            help='number of parallel processes')
    parser_run.set_defaults(func=run_validation)

    # parser for the "report" command
    parser_report = subparsers.add_parser('report', help='generate figures from generated data')
    parser_report.add_argument('--file',
                               action='store',
                               dest='reportfile',
                               default='',
                               help='file name of generated report')
    parser_report.set_defaults(func=report_validation)

    # parse and act
    logging.basicConfig(level=logging.INFO)
    args = parser.parse_args()
    args.func(args)
