import unittest
import os


test_modules = [
	'tests.' + f[:-len('.py')]
	for f in os.listdir(os.path.dirname(__file__))
	if f.startswith('test_') and f.endswith('.py')
]


suite = unittest.TestSuite()
for mod in test_modules :
	suite.addTest(unittest.defaultTestLoader.loadTestsFromName(mod))

unittest.TextTestRunner().run(suite)
