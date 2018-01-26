import unittest


def AbstractTestCase(name, cls) :
	"""Support tests for abstrcat base classes."""
	class BaseTestCase(unittest.TestCase) :
		def run(self, *args, **opts) :
			"""Run the test case only for non-abstract subclasses."""
			if getattr(self, name).__abstractmethods__ :
				return
			else :
				return super(BaseTestCase, self).run(*args, **opts)
	setattr(BaseTestCase, name, cls)
	return BaseTestCase
