"""Event tests

These tests specify the behavior of the Event class.
"""
import unittest
import stocal


class TestEvent(unittest.TestCase) :
	def test_singular_event(self) :
		"""fires at given time for any earlier time and inf after"""
		t = 0.1
		event = stocal.Event({}, {'a':1}, t)
		self.assertEqual(event.next_occurrence(0.), t)
		self.assertEqual(event.next_occurrence(2*t), float('inf'))

	def test_periodic_event(self) :
		"""fires repeatedly with given offset and period"""
		offset = 0.1
		dt = 1.
		event = stocal.Event({}, {'a':1}, offset, dt)
		for t in xrange(100) :
			self.assertAlmostEqual(event.next_occurrence(t), t+offset)


if __name__ == '__main__' :
	unittest.main()
