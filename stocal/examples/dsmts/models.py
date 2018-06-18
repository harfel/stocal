#! /usr/bin/env python3
"""Perform statistical validation tests of the stocal library
"""
import abc
import stocal
try:
    import numpy as np
except ImportError:
    logging.error("dsmts.models suite requires numpy.")
    sys.exit(1)


from stocal._utils import with_metaclass


class DSMTS_Test(with_metaclass(abc.ABCMeta, object)):
    tmax = 50.
    dt = 1.

    def __call__(self, sampler, dt=0, tmax=None):
        """Sample species along sampler every dt time units

        species is a list of species labels that should be sampled.

        Returns a tuple of two elements, the first is a list of all firing
        times, the second a dictionary that holds for each species the
        list of copy numbers at each corresponding time point. If dt is
        given, it specifies the interval at which the sampler is sampled.
        """
        def every(sampler, dt, tmax):
            sampler.tmax = sampler.time
            while sampler.time < tmax:
                transitions = {}
                if sampler.steps and sampler.step >= sampler.steps:
                    break
                sampler.tmax += dt
                for trans in sampler:
                    transitions[trans] = transitions.get(trans, 0) + 1
                yield transitions

        dt = dt or self.dt
        tmax = tmax if tmax is not None else self.tmax

        times = [sampler.time]
        counts = {s: np.array([sampler.state[s]]) for s in self.species}
        it = every(sampler, dt, tmax) if dt else iter(sampler)
        for _ in it:
            times.append(sampler.time)
            for s in self.species:
                counts[s] = np.append(counts[s], [sampler.state[s]])
        return np.array(times), counts

    @abc.abstractproperty
    def process(self): pass

    @abc.abstractproperty
    def initial_state(self): pass

    @abc.abstractproperty
    def species(self): pass

    @classmethod
    def reported_means(cls):
        import os
        from csv import DictReader
    
        dirname = os.path.dirname(__file__)
        fname = cls.__name__ + '-mean.csv'
        path = os.path.join(dirname, fname)
        reader = DictReader(open(path))
        species = reader.fieldnames[1:]
        times = []
        mean = {s: [] for s in species}
        for record in reader:
            times.append(record['time'])
            for s in species:
                mean[s] = np.append(mean[s], [float(record[s])])
        return np.array(times), mean

    @classmethod
    def reported_stdevs(cls):
        import os
        from csv import DictReader
        try:
            import numpy as np
        except ImportError:
                logging.error("DSMTS_Test.reported_stdevs requires numpy.")
                sys.exit(1)

        dirname = os.path.dirname(__file__)
        fname = cls.__name__ + '-sd.csv'
        path = os.path.join(dirname, fname)
        reader = DictReader(open(path))
        species = reader.fieldnames[1:]
        times = []
        mean = {s: [] for s in species}
        for record in reader:
            times.append(record['time'])
            for s in species:
                mean[s] = np.append(mean[s], [float(record[s])**.5])
        return np.array(times), mean


class DSMTS_001_01(DSMTS_Test):
    species = ['X']
    process = stocal.Process([
        stocal.MassAction(['X'], ['X', 'X'], 0.1),
        stocal.MassAction(['X'], [], 0.11)])
    initial_state = {'X': 100}


class DSMTS_001_03(DSMTS_001_01):
    process = stocal.Process([
        stocal.MassAction(['X'], ['X', 'X'], 1.),
        stocal.MassAction(['X'], [], 1.1)])


class DSMTS_001_04(DSMTS_001_01):
    initial_state = {'X': 10}


class DSMTS_001_05(DSMTS_001_01):
    initial_state = {'X': 10000}


class DSMTS_001_07(DSMTS_001_01):
    species = ['X', 'Sink']
    process = stocal.Process([
        stocal.MassAction(['X'], ['X', 'X'], 0.1),
        stocal.MassAction(['X'], ['Sink'], 0.11)])


class DSMTS_002_01(DSMTS_Test):
    species = ['X']
    process = stocal.Process([
        stocal.MassAction([], ['X'], 1.),
        stocal.MassAction(['X'], [], 0.1)])
    initial_state = {}


class DSMTS_002_02(DSMTS_002_01):
    process = stocal.Process([
        stocal.MassAction([], ['X'], 10.),
        stocal.MassAction(['X'], [], 0.1)])


class DSMTS_002_03(DSMTS_002_01):
    # The original test tests for the overloading of global parameters
    # by local parameters. stocal does not have these concepts.
    # We merely check whether the process produces reported results
    # for an immigration rate of 5.
    process = stocal.Process([
        stocal.MassAction([], ['X'], 5.),
        stocal.MassAction(['X'], [], 0.1)])


class DSMTS_002_04(DSMTS_002_01):
    process = stocal.Process([
        stocal.MassAction([], ['X'], 1000.),
        stocal.MassAction(['X'], [], 0.1)])


class DSMTS_002_06(DSMTS_002_01):
    species = ['X', 'Sink']
    process = stocal.Process([
        stocal.MassAction([], ['X'], 10.),
        stocal.MassAction(['X'], ['Sink'], 0.1)])


class DSMTS_002_09(DSMTS_002_01):
    """DSMTS-002-09

    stocal has no equivalent to SBML Events (stocal.transitions.Event
    is something different). Instead, the model's run method samples
    times in two blocks, resetting the 'X' between the two blocks.
    """
    def __call__(self, sampler, dt=0, tmax=None):
        # sample until t=24.
        times0, counts0 = super(DSMTS_002_09, self).__call__(sampler, tmax=24.)

        # advance sampler to t=25.
        sampler.tmax = 25.
        for _ in sampler:
            pass
        
        # update state
        sampler.update_state({'X': 50})

        # sample until t=26.
        times1, counts1 = super(DSMTS_002_09, self).__call__(sampler, tmax=50.)

        # merge state count dictionaries
        counts = {
            s: np.append(counts0[s], counts1[s])
            for s in set(counts0).union(counts1)
        }
        return np.append(times0, times1), counts



class DSMTS_002_10(DSMTS_002_01):
    """DSMTS-002-10

    stocal has no equivalent to SBML Events (stocal.transitions.Event
    is something different). Instead, the model's run method samples
    times in two blocks, resetting the 'X' between the two blocks.
    """
    def __call__(self, sampler, dt=0, tmax=None):
        # sample until t=22.
        times0, counts0 = super(DSMTS_002_10, self).__call__(sampler, tmax=22.)

        # advance sampler to t=22.5
        sampler.tmax = 22.5
        for _ in sampler:
            pass
        
        # update state
        sampler.update_state({'X': 20})

        # advance sampler to t=23.
        sampler.tmax = 23.
        for _ in sampler:
            pass

        # sample until t=26.
        times1, counts1 = super(DSMTS_002_10, self).__call__(sampler, tmax=50.)

        # merge state count dictionaries
        counts = {
            s: np.append(counts0[s], counts1[s])
            for s in set(counts0).union(counts1)
        }
        return np.append(times0, times1), counts


class DSMTS_003_01(DSMTS_Test):
    species = ['P', 'P2']
    process = stocal.Process([
        stocal.MassAction(['P', 'P'], ['P2'], 0.001),
        stocal.MassAction(['P2'], ['P', 'P'], 0.01)])
    initial_state = {'P': 100}


class DSMTS_003_02(DSMTS_003_01):
    process = stocal.Process([
        stocal.MassAction(['P', 'P'], ['P2'], 0.0002),
        stocal.MassAction(['P2'], ['P', 'P'], 0.004)])
    initial_state = {'P': 1000}


class DSMTS_004_01(DSMTS_Test):
    species = ['X']
    process = stocal.Process([
        stocal.MassAction([], {'X': 5}, 1.),
        stocal.MassAction(['X'], [], 0.2)])
    initial_state = {}


class DSMTS_004_02(DSMTS_004_01):
    process = stocal.Process([
        stocal.MassAction([], {'X': 10}, 1.),
        stocal.MassAction(['X'], [], 0.4)])


class DSMTS_004_03(DSMTS_004_01):
    process = stocal.Process([
        stocal.MassAction([], {'X': 100}, 1.),
        stocal.MassAction(['X'], [], 4.)])
