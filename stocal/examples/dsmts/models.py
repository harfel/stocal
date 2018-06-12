#! /usr/bin/env python3
"""Perform statistical validation tests of the stocal library
"""
import abc
import stocal

from stocal._utils import with_metaclass


class DSMTS_Test(with_metaclass(abc.ABCMeta, object)):
    def __call__(self, Algorithm):
        import numpy as np # XXX
        # XXX from validation.run_simulation
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
    
        trajectory = Algorithm(self.process, self.initial_state, tmax=50)
        return sample(trajectory, self.species, dt=1)

        # XXX old __call__ method -- unused
        sampler = Algorithm(self.process, self.initial_state, tmax=0)

        data = []
        while sampler.time < 50:
            sampler.tmax += 1
            for trans in sampler: pass
            data.append([sampler.state[s] for s in self.species])
        return data

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
        try:
            import numpy as np
        except ImportError:
                logging.error("DSMTS_Test.reported_means requires numpy.")
                sys.exit(1)
    
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
    # instead of using Event's this simply splits sampling into two
    # stages, resetting species 'X' after the first stage.
    def XXX___call__(self, Algorithm):
        sampler = Algorithm(self.process, self.initial_state, tmax=0)

        data = []
        while sampler.time < 25:
            sampler.tmax += 1
            for trans in sampler: pass
            data.append([sampler.state[s] for s in self.species])
        sampler.update_state(X=50)
        while sampler.time < 50:
            sampler.tmax += 1
            for trans in sampler: pass
            data.append([sampler.state[s] for s in self.species])
        return data


# XXX 002_10


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

# XXX 003_03
# XXX 003_04
# XXX 003_05 ?


class DSMTS_004_01(DSMTS_Test):
    species = ['X']
    process = stocal.Process([
        stocal.MassAction([], {'X': 5}, 1.),
        stocal.MassAction(['X'], [], 0.2)])
    initial_state = {}


class DSMTS_004_02(DSMTS_004_01):
    process = stocal.Process([
        stocal.MassAction([], {'X': 10}, 1.),
        stocal.MassAction(['X'], [], 0.2)])


class DSMTS_004_03(DSMTS_004_01):
    process = stocal.Process([
        stocal.MassAction([], {'X': 100}, 1.),
        stocal.MassAction(['X'], [], 0.2)])
