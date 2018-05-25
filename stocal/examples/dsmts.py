#! /usr/bin/env python3
"""Perform statistical validation tests of the stocal library
"""
import abc
import logging
import os
import stocal

from stocal._utils import with_metaclass


logging.basicConfig(level=logging.INFO)

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


class DSMTS_001_01(DSMTS_Test):
    species = ['X']
    process = stocal.Process([
        stocal.MassAction(['X'], ['X', 'X'], 0.1),
        stocal.MassAction(['X'], [], 0.11)])
    initial_state = {'X': 100}


def accumulate(sample, mean=[[0.] for i in range(50)], runs=[0]):
    # XXX do this in stocal.examples.validation
    for i,data in enumerate(mean):
        for j,species in enumerate(sample[i]):
            mean[i][j]=(runs[0]*data[j]+species)/(runs[0]+1)
    runs[0] += 1
    return mean


def main():
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


if __name__ == '__main__':
    main()
