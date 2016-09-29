# FATSLiM Changelog
All notable changes to FATSLiM are documented in this file.
Version numbers comply with [Python recommandation](https://www.python.org/dev/peps/pep-0440/).

## [0.2.1] - Unreleased
### Added
- Description of the "--idfreq" option in the documentation

### Changed
- Default frequence for membrane identification is now 1 (every frame)

### Fixed
- Bug in lipid direction calculation
- Lack of verbose output for "membranes" command
- Ignored "--idfreq" option leading to membrane identified only once per trajectory
- Bad test file for thickness command xvg export
- --end-frame option interpreted as float an raising error when used

## [0.2.0] - 2016-07-27
### Added
- Full documentation.
- Updated thickness calculation algorithm
- --begin-frame/--end-frame options to select begining/ending frame index

### Changed
- --begin and --end options new refers to timesteps (same as GROMACS) rather than frame index

### Fixed
- Missing "--interacting-group"

## [0.1.2] - 2016-03-21
### Added
- Changelog file.
- [Travis CI](https://travis-ci.org/FATSLiM/fatslim) support.
- [Coveralls](https://coveralls.io/github/FATSLiM/fatslim) support

### Fixed
- README.rst is now correctly rendered by github.
- Bug in leaflet attribution routine.


## [0.1] - 2016-02-01
First public release

[Unreleased]: https://github.com/FATSLiM/fatslim/tree/develop
[0.1]: https://github.com/FATSLiM/fatslim/releases/tag/v0.1