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
                mean[s].append(float(record[s]))
        return times, mean

    @classmethod
    def reported_stdevs(cls):
        import os
        from csv import DictReader
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
                mean[s].append(float(record[s])**.5)
        return times, mean


class DSMTS_001_01(DSMTS_Test):
    species = ['X']
    process = stocal.Process([
        stocal.MassAction(['X'], ['X', 'X'], 0.1),
        stocal.MassAction(['X'], [], 0.11)])
    initial_state = {'X': 100}


def main():
    import logging
    import os

    logging.basicConfig(level=logging.INFO)

    # XXX do this in stocal.examples.validation (if at all)
    def download_dsmts(url='https://github.com/darrenjw/dsmts/archive/master.zip'):
        from urllib.request import urlopen
        from zipfile import ZipFile
        from io import BytesIO
        logging.info("Downloading dsmts")
        ZipFile(BytesIO(urlopen(url).read())).extractall()

    dsmts_downloaded = False

    if not os.path.exists('dsmts-master'):
        download_dsmts()
        dsmts_downloaded = True

    # collect data for each algorithm and model
    for n in range(10):
        sample = DSMTS_001_01()(stocal.algorithms.DirectMethod)
        mean = accumulate(sample)
        print('\n'.join(('%d '%t) + ' '.join(str(s) for s in row) for t, row in enumerate(mean, 1)))
        print()
    # perform statistical tests
    # generate report

    if dsmts_downloaded:
        pass # delete dsmts-master directory

if __name__ == '__main__' :
    m = DSMTS_001_01()
    print m.reported_means()
