import unittest
import os

os.chdir(os.path.dirname(__file__))

test_modules = [
	f[:-len('.py')]
	for f in os.listdir('.')
	if f.startswith('test_') and f.endswith('.py')
]


suite = unittest.TestSuite()
for mod in test_modules :
	suite.addTest(unittest.defaultTestLoader.loadTestsFromName(mod))

unittest.TextTestRunner().run(suite)
