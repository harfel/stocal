# Changelog

## [2.0]

### Removed
- Removed stocal.algorithm.TrajectorySampler
- Removed stocal.algorithms.AndersonNRM
- Removed stocal.Process.trajectory
- Removed stocal.ReactionRule
- Removed keyword argument t to StochasticSimulationAlgorithm
- Removed stocal.experimental.samplers

### Changed
- TransitionRule.infer_transitions can no longer be called with dicts, use multisets instead
- PriorityQueue.keys and MultiDict.keys return generators rather than lists



## [1.2] - 2018-08-10

### Added
- Added Gibson & Bruck's NextReactionMethod
- Added Cao et al's tau leaping method in stocal.experimental.tauleap
- Added statistical validation suite in stocal.examples.validation and stocal.samples.dsmts
- Added modular trajectory sampling interface stocal.experimental.samplers
- Added flattening of rule-based processes into static processes
- StochasticSimulationAlgorithm instances accept an optional random seed

### Deprecated
- Deprecated AndersonNRM in favour of AndersonMethod
- Deprecated TrajectorySampler in favour of StochasticSimulationAlgorithm
- Deprecated ReactionRule in favour of TransitionRule
- Deprecated Process.trajectory in favour of Process.sample

### Changed
- TrajectorySampler instances now use an internal random number generator (sampler.rng) rather than python's global one.
- Improved performance of DirectMethod and AndersonMethod
- Improved performance of transition inference in TransitionRule

## [1.1.2] - 2018-05-23

### Fixed
- Fixed issue #4 in FirstReactionMethod

### Changed
- stocal.examples.type_rules now does what it claims to do
- Events that occurr at tmax are now explicitly included in a trajectory


## [1.1.1] - 2018-03-21

### Fixed
- Fixed issue 3 in AndersonNRM

### Changed
- MassAction.__repr__ now also prints rate constant


## [1.1] - 2018-03-08

### Fixed
- Fixed issue 1 in Event.next_occurrence
- Fixed issue 2 in ReactionRule.infer_transitions
- Fixed issue with printing transitions with non-string reactants

### Added
- Added type support of reaction rules via ReactionRule.signature
- Added support for time dependent reaction rates

### Changed
- ReactionRule.Transition can alternatively be specified through return type annotation in python3


## [1.0.2] - 2018-02-07

### Added
- Added python3 compatibility

### Deprecated
- Deprecated dict support of several functions in favor of multisets
- Deprecated some method parameter names

### Removed
- Removed requirement for ReactionRule.order

### Changed
- Changed package organization


## [1.0.1] - 2018-01-29

### Fixed
- Fixed error in tutorial


## [1.0] - 2018-01-24

initial release
