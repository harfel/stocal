#! /usr/bin/env python

from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()


setup(name = "stocal",
      version = "1.0.2",
      description = "simple rule-based stochastic simulation",
      long_description = readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
      ],
      url = "https://github.com/harfel/",
      author = "Harold Fellermann",
      author_email = "harold.fellermann@newcastle.ac.uk",
      license='MIT',
      packages = ["stocal", "stocal.examples", "stocal.tests"],
      include_package_data=True,
      zip_safe = True,
      test_suite = 'stocal.tests')
