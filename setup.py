#! /usr/bin/env python

from distutils.core import setup

setup(
	name = "stocal",
	version = "1.0.1",
	description = "simple rule-based stochastic simulation",
	author = "Harold Fellermann",
	author_email = "harold.fellermann@newcastle.ac.uk",
	url = "https://github.com/harfel/",
	packages = ["stocal", "stocal.examples"],
	package_data = {'stocal' : ['doc/*.html']},
)
