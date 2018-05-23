# Changelog

## [1.2]

### Added
- Added modular trajectory sampling interface stocal.experimental.samplers
- TrajectorySampler instances accept an optional random seed
- Added Gibson & Bruck's NextReactionMethod

### Changed
- TrajectorySampler instances now use an internal random number generator (sampler.rng) rather than python's global one.
- Improved performance of DirectMethod and AndersonNRM


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
