#! /usr/bin/env python
import os
from setuptools import setup

here = os.path.abspath(os.path.dirname(__file__))

def find_version(*file_paths):
    import re
    def read(*parts):
        import codecs
        with codecs.open(os.path.join(here, *parts), 'r') as fp:
            return fp.read()

    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


def readme():
    with open('README.md') as f:
        return f.read()


setup(name = "stocal",
      version = find_version("stocal", "__init__.py"),
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
      packages = ["stocal", "stocal.examples", "stocal.experimental", "stocal.tests"],
      include_package_data=True,
      zip_safe = True,
      test_suite = 'stocal.tests',
      install_requires=['pqdict'])
