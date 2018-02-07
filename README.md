# stocal

A python framework for regular and rule-based stochastic simulations.


## What is Stocal?

Stocal is a framework for stochastic simulation of continuous
time Markov processes (also known as [Gillespie simulations](https://en.wikipedia.org/wiki/Gillespie_algorithm)). 

Features of stocal:
* support for reactions of any order
* support for unique and periodic (deterministic) events
* support for rules that generate novel reactions on the fly
* support for complex chemical states

Stocal works with python version 2.7 as well as version 3.5.


### Basic Usage

Running a simple stochastic simulation is straight forward:
```python
import stocal

# Define a stochastic process
process = stocal.Process([
	stocal.MassAction({'A': 2}, {'A2': 1}, 1.0),
	stocal.MassAction({'A2': 1}, {'A': 2}, 10.0),
	stocal.Event({}, {'A': 100}, 0.0, frequency=10.0),
])

# Sample a stochastic trajectory of the process
initial_state = {}
trajectory = process.trajectory(initial_state, tmax=100) :
	print trajectory.time, trajectory.state['A2']
```


### API / Tutorial

A tutorial on how to use stocal can be found [here](https://github.com/harfel/stocal/wiki).
Various usage examples are provided in stocal/examples.
The package API is thoroughly documented and can be accessed through
pydoc. The behavior of stocal is specified via tests. The test suite
can be run with `python setup.py test`.


### Installation

The latest stable release of stocal is available from the python package index:
```shell
pip install stocal
```
The development version can be obtained from github using the following commands:
```shell
git clone https://github.com/harfel/stocal.git
cd stocal
git checkout develop
sudo python setup.py install
```


## Issue Reporting and Contributing

Please post any issues that might occur with stocal on the
<a href="https://github.com/harfel/stocal/issues">github issue
tracker</a>.
If you are interested in contributing to stocal, pull requests
and any other inquiries will be dealt with as soon as possible.


## License

Stocal is distributed under the MIT Software license.
(c) 2018 Harold Fellermann
