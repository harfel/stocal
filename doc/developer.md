# stocal developer guide

The following are guidelines for developers of stocal. If you are
interested in using stocal as a library/application, read the tutorial
instead.


## Contributing

Pull-requests are always welcome and will be considered as soon as
possible. Your request is more likely to be accepted if it adheres to
the guidelines detailed in this document.


## Project development

stocal uses git for distributed version control and follows the gitflow
workflow. Pull requests should only be made for feature, bugfix, hotfix,
or support branches. Pull requests to the master or develop branch are
likely to be rejected straight away.


## Style
If in doubt, pylint decides.


## Tests
The test suite for stocal is contained in `stocal.tests` and can be
run using
```bash
python -m unittest discover stocal.tests
python3 -m unittest discover stocal.tests
```

Passing of all tests is an enforced requirement for all code merged
into the develop and master branches.

Coverage analysis can be performed using the
[third party](https://pypi.python.org/pypi/coverage) `coverage` module.
The command lines are:
```bash
coverage run --source=stocal --omit='stocal/tests/*' setup.py test
coverage html
```

New functionality should be accompanied by tests. For novel
implementations of defined interfaces (abstract classes), such as
Transition or Rule, stocal.test.abstract_test offers an infrastructure
to derive implementation test cases from interface test cases. See
`pydoc stocal.tests` for more information.


## Validation
Stocal ships with a validation suite for stochastic simulation
algorithms. The validation suite is based on the discrete stochastic
simulation model test suite DSMTS. To run validations, call
$ python stocal/examples/validation.py run N
from the command line. To generate a validation report, run
$ python stocal/examples/validation.py report
This generates a file validation.tex that can be compiled with pdflatex.
See
$ python stocal/examples/validation.py -h
for more information. The DSMTS user guide recommends N=1,000 as an
absolute minimum to be run regularly in conjunction with unit test,
and n=100,000 or n=1,000,000 for a thorough statistical analysis.
A rudimentary LaTeX template that collates report results into a
single document can be found in doc/validation.tex


## Releases

When preparing a new release, these steps should be followed

 * git flow release start
 * ensure an optimal code coverage of the test suite
 * ensure that all tests pass
 * ensure that any novel algorithm passes validation
 * ensure that documentation (README, tutorial, etc.) is up to date
 * update CHANGELOG.md
 * bump the version number
 * git flow release finish
 * sudo setup.py sdist upload
 * publish release on github
