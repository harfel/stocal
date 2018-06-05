#! /usr/bin/env python3
"""Perform statistical validation tests of the stocal library
"""
import abc
import stocal

from stocal._utils import with_metaclass


class DSMTS_Test(with_metaclass(abc.ABCMeta, object)):
    def __call__(self, Algorithm):
        sampler = Algorithm(self.process, self.initial_state, tmax=0)

        data = []
        while sampler.time < 50:
            sampler.tmax += 1
            for trans in sampler:
                pass
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
