# Copyright 2018-2020 Harold Fellermann
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""Facilities to work with stochastic trajectories"""

class Trajectory(object):
    """Access to a single random sample of a stochastic process

    Trajectory instances contain a time access (trajectory.times) at which the
    system state is sampled, and a dictionary of copy number axes
    (trajectory.species), one for each species that occurred during simulation.

    Rather than directly instantiating a Trajectory, user code should call
    Process.trajectory() to obtain an instance.
    """
    def __init__(self, sampler, times=None):
        """XXX init doc"""
        self.sampler = sampler
        self.species = {}
        self.times = times
        # XXX make trajectory lazy
        self._sample()
    
    def _sample(self):
        if self.times is not None:
            for idx, t in enumerate(self.times):
                self.sampler.until(t)
                for s in self.sampler.state:
                    if s not in self.species:
                        self.species[s] = len(self.times)*[0]
                    self.species[s][idx] = self.sampler.state[s]

        else:
            self.times = []
            for trans in self.sampler:
                for s in self.sampler.state:
                    if s not in self.species:
                        self.species[s] = len(self.times)*[0]
                    self.species[s].append(self.sampler.state[s])
                for s in self.species:
                    if s in self.sampler.state:
                        continue
                    self.species[s].append(0)
                self.times.append(self.sampler.time)



if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    from stocal import *

    r1 = MassAction({'A': 2}, {'A2': 1}, 1.0)
    r2 = MassAction({'A2': 1}, {'A': 2}, 10.0)
    process = Process([r1, r2])
    
    trajectory = process.trajectory({'A': 100}, tmax=1.)
    for species in trajectory.species:
        plt.step(trajectory.times, trajectory[species])
    #plt.show()

    times = np.linspace(0., 1., 11)
    trajectory = process.trajectory({'A': 100}, times=times)
    for species in trajectory.species:
        plt.step(times, trajectory[species])
    #plt.show()
