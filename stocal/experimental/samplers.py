# Copyright 2018 Harold Fellermann
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
"""Chainable process samplers

This module is currently experimental. Its API might change without
notification, and there is no guarantee that it will be adopted by
the stable stocal core.

The samplers module provides a convenient framework to sample
stocal.Process instances. The framework is based on Sampler instances
that can be chained to customize the way a trajectory is sampled. The
simplest way to obtain a Sampler for a given process is to call the
processes sample method:

>>> sampler = process.sample(initial_state, tstart=0)

Process.sample is a factory method that selects and initializes an
appropriate stochastic simulation algorithm depending on the properties
of the given process.

The sampler obtained this way simply gives access to each iteration of
the algorithm: 

>>> for time, state, transitions in sampler:
>>>     print(time, state)

The yielded triple gives information of system time and state, plus
a dictionary that counts all transitions that occurred since the last
iteration (having only one key with value one in the above scenario).


All samplers expose methods that fine-tune the way a trajectory is
sampled. For example, to sample a trajectory until a given stop time
is reached, use:

>>> process.sample(initial_state).until(time=stop_time)

To sample only the first n steps, use:

>>> process.sample(initial_state).until(steps=n)

This will return samplers augmented with the given stop criteria.


By default, samplers yield one datapoint for each Transition generated
by the simulation algorithm. To only yield datapoints after m steps,
use:

>>> process.sample(initial_state).every(steps=m).until(steps=n)

or, to return a data point every dt time units, use:

>>> process.sample(initial_state).every(time=dt).until(steps=n)


Iterating over Sampler.every will yield the system state at the given
intervals. If you instead want to obtain the average system state over
the sampled interval, iterate over:

>>> process.sample(initial_state).average(steps=m)

or

>>> process.sample(initial_state).average(time=dt)

Now, the yielded state is the average state over the entire sampled
interval. In the first case, each state has the same weight in the
average. In the second case, states are weighed by the time they
occupied in the sampled interval.


To obtain data points only after certain Transition's have taken place,
use:

>>> process.sample(initial_state).filter([tran_a, trans_b, ...])


Note that the order in which Sampler's are chained is significant:

>>> process.sample(initial_state).average(steps=10).filter([trans])

will average the system states in blocks of ten and yield only those
blocks that feature trans.

>>> process.sample(initial_state).filter([trans]).average(steps=10)

in contrast, will take only system states ater each trans transition
and calculate averages of these states.

See the Sampler documentation for a full explanation of Sampler factory
methods.

# XXX Chaining of Samplers as described in the doc is currently broken
"""
import abc
try:
    from itertools import izip as zip
    range = xrange
except ImportError:
    pass

import stocal.transitions
from stocal._utils import with_metaclass
from stocal.structures import multiset

StandardProcess = stocal.transitions.Process

class Process(StandardProcess):
    """Process class that supports the new sampler protocol
    
    Importing the module performs a monkey patch that replaces
    stocal.transitions.Process with the samplers.Process class.
    """
    def sample(self, state, tstart=0., tmax=float('inf'), steps=None, every=None, seed=None):
        """Create trajectory sampler for given state
        
        The method selects a suitable stochastic simulation algorithm
        given the Process'es Transition's and Rule's, and wraps this
        algorithm into a Sampler, specified via optional arguments.
        If tmax is supplied, the sampler will terminate when the sampler
        reaches the given time. If steps is specified, the sampler
        terminates after the given number of steps instead. Both
        stop criteria can be supplied to stop at whichever event occurs
        first. If tmax or steps (bot not both) are given, every can be
        used to specify a sampling interval. In conjunction with tmax,
        the sampler will yield trajectory states in given time intervals.
        In conjunction with steps, the sampler will yield results after
        a given number of steps.
        """
        algorithm = super(Process, self).trajectory(state, tstart=tstart, seed=seed)
        sampler = _Wrapper(algorithm)

        # instantiate requested sampling method
        if every is not None:
            if tmax != float('inf') and steps is not None:
                raise ValueError("every can only be provided when either steps or tmax is given, not both.")
            elif tmax != float('inf'):
                return sampler.every(time=every).until(time=tmax)
            elif steps is not None:
                if not isinstance(every, int):
                    raise ValueError("every must be an integer in combination with steps.")
                return sampler.every(steps=every).until(steps=steps)
            else:
                raise ValueError("every can only be used in conjunction with a stop criterion.")
        else:
            if tmax != float('inf') and steps is not None:
                return sampler.until(time=tmax, steps=steps)
            elif tmax != float('inf'):
                return sampler.until(time=tmax)
            elif steps is not None:
                return sampler.until(steps=steps)
            else:
                return sampler


# monkey patch stocal.Process
stocal.Process = Process
stocal.transitions.Process = Process


class Sampler(with_metaclass(abc.ABCMeta, object)):
    """Abstract base class for Samplers
    
    Concretizations have to provide an __iter__ method.
    """
    def __init__(self, sampler):
        self.sampler = sampler
        self.algorithm = sampler.algorithm

    @abc.abstractmethod
    def __iter__(self):
        """Iterator support.
        
        All Sampler classes need to provide iterators that yield
        information about the sampled trajectory, namely a triple
        (time, state, transitions), where time and state are the
        system time and state and transitions is a dictionary counting
        all transitions having occurred within one iteration.
        """
        raise StopIteration

    def __getattr__(self, attr):
        """Delegate attribute lookup through sampler chain."""
        return getattr(self.sampler, attr)

    def until(self, time=float('inf'), steps=None):
        """Return sampler with given stop criterion."""
        if time != float('inf') and steps != None:
            return UntilTimeSampler(self, time).until(steps=steps)
        elif time != float('inf'):
            return UntilTimeSampler(self, time)
        elif steps:
            return UntilStepSampler(self, steps)
        else:
            raise ValueError("Either time or steps must be given.")

    def every(self, time=float('inf'), steps=None):
        """Return sampler that yields every given time or steps."""
        if time != float('inf') and steps != None:
            raise ValueError("time and steps cannot both be given.")
        elif time != float('inf'):
            return EveryTimeSampler(self, time)
        elif steps:
            return EveryStepSampler(self, steps)
        else:
            raise ValueError("Either time or steps must be given.")

    def average(self, time=float('inf'), steps=None):
        """Return sampler that averages over given time or steps."""
        if time != float('inf') and steps != None:
            raise ValueError("time and steps cannot both be given.")
        elif time != float('inf'):
            return AverageTimeSampler(self, time)
        elif steps:
            return AverageStepSampler(self, steps)
        else:
            raise ValueError("Either time or steps must be given.")

    def filter(self, transitions):
        """Return sampler that yields only if any given transition occurred."""
        return FilteredSampler(self, transitions)


class _Wrapper(Sampler):
    """Wrapper to use StochasticSimulationAlgorithm with Sampler interface

    This class only exists for transition to the new interface and will
    be removed when no longer necessary. Please do not use explicitly
    in productive code.
    """
    def __init__(self, sampler):
        self.sampler = sampler
        self.algorithm = sampler

    def __iter__(self):
        traj = self.sampler
        for transition in traj:
            yield traj.time, traj.state, { transition: 1}


class UntilTimeSampler(Sampler):
    """Sample until time equals given final time

    Iterate over the underlying sampler until sampler.time
    equals the given final time.
    """
    def __init__(self, sampler, time):
        super(UntilTimeSampler, self).__init__(sampler)
        self.tmax = time

    def __iter__(self):
        while True:
            time, transition, args = self.propose_potential_transition() # TODO: needs to iterate over self.sampler!

            if time > self.tmax:
                break
            else:
                self.perform_transition(time, transition, *args)
                yield time, self.state, transition
        self.algorithm.time = self.tmax


class UntilStepSampler(Sampler):
    """Sample for given number of steps

    Iterate over the underlying sampler for a given number of steps.
    """
    def __init__(self, sampler, steps):
        super(UntilStepSampler, self).__init__(sampler)
        self.steps = steps

    def __iter__(self):
        for n, data in zip(range(self.steps), self.sampler):
            yield data


class EveryTimeSampler(Sampler):
    """Yield samples every dt time units

    Iterates over the underlying sampler and returns time, state
    and transitions every dt time units. If initialized with skip=True,
    the sampler skips time points where no system change has occurred.
    """
    def __init__(self, sampler, time, skip=False):
        super(EveryTimeSampler, self).__init__(sampler)
        self.dt = time
        self.skip = skip

    def __iter__(self):
        algorithm = self.algorithm
        time = self.sampler.time
        transitions = multiset()
        while True:
            ptime, trans, args = algorithm.propose_potential_transition() # TODO: needs to iterate over self.sampler!
            if ptime > time+self.dt:
                time += self.dt
                yield time, self.state, transitions
                transitions = multiset()

                if ptime == float('inf'):
                    break
                while self.skip and time+self.dt < ptime:
                    time += self.dt

            transitions[trans] += 1
            self.algorithm.perform_transition(ptime, trans, *args)


class EveryStepSampler(Sampler):
    """Yield samples every n steps

    Iterates over the underlying sampler and returns time, state
    and transitions every n steps.
    """
    def __init__(self, sampler, steps):
        super(EveryStepSampler, self).__init__(sampler)
        self.steps = steps

    def __iter__(self):
        while True:
            transitions = multiset()
            for n, data in zip(range(self.steps), self.sampler):
                transitions += data[2]
            if not transitions:
                break
            yield self.time, self.state, transitions


class AverageTimeSampler(EveryTimeSampler):
    """Average samples over dt time units

    Iterates over the underlying sampler and averages the system state
    over the given period of time. When averaging, each state is weighed
    by the time that it persisted.
    """
    def __iter__(self):
        algorithm = self.algorithm
        time = self.sampler.time
        transitions = multiset()
        averages = multiset()
        while True:
            ptime, trans, args = algorithm.propose_potential_transition() # TODO: needs to iterate over self.sampler!

            if ptime > time+self.dt:
                time += self.dt
                yield time, 1./self.dt*averages, transitions
                if ptime == float('inf'):
                    break
                transitions = multiset()
                averages = multiset()
                while ptime > time+self.dt:
                    time += self.dt
                    if not self.skip:
                        yield time, self.state, {}

            transitions[trans] += 1
            averages += (ptime-max(time, algorithm.time))*algorithm.state
            algorithm.perform_transition(ptime, trans, *args)


class AverageStepSampler(EveryStepSampler):
    """Average samples over dt time units

    Iterates over the underlying sampler and averages system state
    over the given number of steps. Each state has equal qeight in the
    average.
    """
    def __iter__(self):
        while True:
            transitions = multiset()
            averages = multiset()
            for n, data in zip(range(self.steps), self.sampler):
                transitions += data[2]
                averages += data[1]
            if not transitions:
                break
            yield self.time, 1./self.steps*averages, transitions


class FilteredSampler(Sampler):
    """Yield samples only when one of the given transitions occurred.

    Iterates over the underlying sampler and returns time, state
    and transitions if one of the given transitions occurred.

    TODO: In addition to Transition instances, this sampler could
    accept Transition classes, Rule instances and Rule classes and
    yield after any transition that equals any given Transition class,
    isinstance of any given Transition class, or has been infered from
    any given Rule instance or class.
    """
    def __init__(self, sampler, transitions):
        super(FilteredSampler, self).__init__(sampler)
        self.transitions = transitions

    def __iter__(self):
        for time, state, transitions in self.sampler:
            if any(trans in transitions for trans in self.transitions):
                yield time, state, transitions
