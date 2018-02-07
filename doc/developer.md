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
python setup.py test
python3 setup.py test
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


## Releases

When preparing a new release, these steps should be followed

 * git flow release start
 * ensure an optimal code coverage of the test suite
 * ensure that all tests pass
 * ensure that documentation (README, tutorial, etc.) is up to date
 * bump the version number
 * git flow release finish
 * sudo setup.py sdist upload
 * publish release on github
