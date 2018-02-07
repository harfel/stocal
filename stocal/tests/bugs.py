import unittest
import stocal

class Issue1(unittest.TestCase):
    """Issue 1 [CLOSED]: https://github.com/harfel/stocal/issues/1

    Opened 07-02-2018
    Closed 07-02-2018
    """
    def test_text_occurrence_delayed(self):
        """assert that repeated events do not fire before first occurrence"""
        event = stocal.Event(['a'], ['z'], 1, 1)
        self.assertEqual(event.next_occurrence(0), 1)


if __name__ == '__main__':
    unittest.main()
