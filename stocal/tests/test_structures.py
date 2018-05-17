"""Tests for the stocal.structures"""
import unittest
from stocal.structures import multiset


class TestMultiset(unittest.TestCase):
    """Tests for structures.multiset"""
    def test_init_assert_number_values(self):
        """raise ValueError for non-numeric values"""
        with self.assertRaises(TypeError):
            multiset({'a': None})
        with self.assertRaises(TypeError):
            multiset({'a': object()})
        with self.assertRaises(TypeError):
            multiset({'a': "string"})

    def test_init_from_dict(self):
        """assert initialization with dictionary"""
        multiset({'a':1, 'z':2})

    def test_init_from_multiset(self):
        """assert initialization from multiset"""
        mset_a = multiset({'a':1, 'z':2})
        mset_b = multiset(mset_a)
        self.assertIsNot(mset_a, mset_b)
        self.assertEqual(mset_a, mset_b)

    def test_init_with_keywords(self):
        """assert initialization with keywords"""
        mset_a = multiset({'a':1, 'z':2})
        mset_b = multiset(a=1, z=2)
        self.assertEqual(mset_a, mset_b)

    def test_init_from_unordered_list(self):
        """assert initialization with unordered list"""
        mset_a = multiset(['z', 'a', 'z'])
        mset_b = multiset({'a':1, 'z':2})
        self.assertEqual(mset_a, mset_b)

    def test_init_from_unordered_list_special_case(self):
        """assert initialization with unordered list ['az']"""
        mset_a = multiset(['az'])
        mset_b = multiset({'az':1})
        self.assertEqual(mset_a, mset_b)

    def test_domain(self):
        """test multiset domain / underlying set"""
        mset = multiset({'a':1, 'z':2})
        self.assertEqual(mset.domain, set(['a', 'z']))

    def test_null(self):
        """assert empty multiset is False"""
        self.assertFalse(multiset())
        self.assertTrue(multiset({'a':1}))

    def test_eq(self):
        """test equality"""
        mset_a = multiset({'a':1, 'z':2})
        mset_b = multiset({'a':1, 'z':2})
        self.assertEqual(mset_a, mset_b)

    def test_eq_zero(self):
        """test equality with and without explicit zero"""
        mset_a = multiset()
        mset_b = multiset({'a':0})
        self.assertEqual(mset_a, mset_b)

    def test_ne(self):
        """test inequality"""
        mset_a = multiset({'a':1, 'z':2})
        mset_b = multiset({'a':1, 'z':1})
        self.assertNotEqual(mset_a, mset_b)

    def test_le(self):
        """test left subset"""
        mset = multiset({'a':1, 'z':2})
        self.assertLessEqual(mset, multiset({'a':1, 'z':2}))
        self.assertLessEqual(mset, multiset({'a':2, 'z':2}))
        self.assertLessEqual(mset, multiset({'a':1, 'x':1, 'z':2}))
        self.assertFalse(mset <= multiset({'a':1, 'z':1}))
        self.assertFalse(mset <= multiset())

    def test_lt(self):
        """test left true subset"""
        mset = multiset({'a':1, 'z':2})
        self.assertFalse(mset < multiset({'a':1, 'z':2}))
        self.assertLess(mset, multiset({'a':2, 'z':2}))
        self.assertLess(mset, multiset({'a':1, 'x':1, 'z':2}))
        self.assertFalse(mset < multiset({'a':1, 'z':1}))
        self.assertFalse(mset < multiset())

    def test_ge(self):
        """test right subset"""
        mset = multiset({'a':1, 'z':2})
        self.assertGreaterEqual(mset, multiset({'a':1, 'z':2}))
        self.assertGreaterEqual(mset, multiset({'a':1, 'z':1}))
        self.assertGreaterEqual(mset, multiset({'a':1}))
        self.assertFalse(mset >= multiset({'a':2, 'z':2}))
        self.assertGreaterEqual(mset, multiset())

    def test_gt(self):
        """test right true subset"""
        mset = multiset({'a':1, 'z':2})
        self.assertFalse(mset > multiset({'a':1, 'z':2}))
        self.assertGreater(mset, multiset({'a':1, 'z':1}))
        self.assertGreater(mset, multiset({'a':1}))
        self.assertFalse(mset >= multiset({'a':2, 'z':2}))
        self.assertGreater(mset, multiset())

    def test_getitem(self):
        """test getitem"""
        mset = multiset({'a':1, 'z':2})
        self.assertEqual(mset['a'], 1)
        self.assertEqual(mset['z'], 2)
        self.assertEqual(mset['x'], 0)

    def test_setitem(self):
        """test setitem"""
        mset = multiset({'a':5})
        mset['a'] = 1
        mset['z'] = 2
        self.assertEqual(mset['a'], 1)
        self.assertEqual(mset['z'], 2)

    def test_setitem_zero(self):
        """test setitem"""
        mset = multiset({'a':5})
        mset['a'] = 0
        self.assertFalse('a' in mset.domain)

    def test_delitem(self):
        """test delitem"""
        mset = multiset({'a':1, 'z':2})
        del mset['z']
        self.assertEqual(mset['z'], 0)

    def test_update(self):
        """test that update deletes elements"""
        mset = multiset({'a':5})
        mset.update(a=0)
        self.assertFalse('a' in mset.domain)

    def test_len(self):
        """test len / number of items"""
        mset = multiset({'a':1, 'z':2})
        self.assertEqual(len(mset), 3)

    def test_contains_element(self):
        """test element in"""
        mset = multiset({'a':1, 'z':2})
        self.assertIn('a', mset)
        self.assertTrue('x' not in mset)

    def test_contains_multiset(self):
        """test subset of"""
        mset = multiset({'a':1, 'z':2})
        self.assertIn(multiset(), mset)
        self.assertIn(mset, multiset({'a':1, 'z':2}))
        self.assertIn(mset, multiset({'a':2, 'z':2}))
        self.assertIn(mset, multiset({'a':1, 'x':1, 'z':2}))

    def test_add(self):
        """test addition"""
        mset_a = multiset({'a':1, 'z':2})
        mset_b = multiset({'a':1, 'x':2})
        mset_c = multiset({'a':2, 'x':2, 'z':2})
        self.assertEqual(mset_a + mset_b, mset_c)
        self.assertNotEqual(mset_a, mset_c)
        self.assertNotEqual(mset_b, mset_c)

    def test_iadd(self):
        """test inplace addition"""
        mset_a = multiset({'a':1, 'z':2})
        mset_b = multiset({'a':1, 'x':2})
        mset_c = multiset({'a':2, 'x':2, 'z':2})
        mset_d = mset_a
        mset_a += mset_b
        self.assertEqual(mset_a, mset_c)
        self.assertIs(mset_a, mset_d)

    def test_union(self):
        """test union"""
        mset_a = multiset({'a':1, 'z':2})
        mset_b = multiset({'a':1})
        mset_c = multiset({'x':2})
        mset_d = multiset({'a':2, 'x':2, 'z':2})
        self.assertEqual(mset_a.union(mset_b, mset_c), mset_d)
        self.assertIsNot(mset_a.union(mset_b, mset_c), mset_a)

    def test_sub(self):
        """test difference"""
        mset_a = multiset({'a':2, 'x':2, 'z':2})
        mset_b = multiset({'a':1, 'z':2})
        mset_c = multiset({'a':1, 'x':2})
        self.assertEqual(mset_a - mset_b, mset_c)
        self.assertNotEqual(mset_a, mset_c)
        self.assertNotEqual(mset_b, mset_c)

    def test_sub_with_overlap(self):
        """test difference"""
        mset_a = multiset({'a':2, 'x':2, 'z':2})
        mset_b = multiset({'a':1, 'b':4, 'z':2})
        mset_c = multiset({'a':1, 'x':2})
        self.assertEqual(mset_a - mset_b, mset_c)
        self.assertNotEqual(mset_a, mset_c)
        self.assertNotEqual(mset_b, mset_c)

    def test_isub(self):
        """test inplace difference"""
        mset_a = multiset({'a':2, 'x':2, 'z':2})
        mset_b = multiset({'a':1, 'x':2})
        mset_c = multiset({'a':1, 'z':2})
        mset_d = mset_a
        mset_a -= mset_b
        self.assertEqual(mset_a, mset_c)
        self.assertIs(mset_a, mset_d)

    def test_difference(self):
        """test difference"""
        mset_a = multiset({'a':2, 'x':2, 'z':2})
        mset_b = multiset({'a':1, 'b':1})
        mset_c = multiset({'x':2})
        mset_d = multiset({'a':1, 'z':2})
        self.assertEqual(mset_a.difference(mset_b, mset_c), mset_d)
        self.assertIsNot(mset_a.difference(mset_b, mset_c), mset_a)

    def test_symmetric_difference(self):
        """test symmetric difference"""
        mset_a = multiset({'a':2, 'x':2, 'z':2})
        mset_b = multiset({'a':1, 'b':1})
        mset_c = multiset({'x':2})
        mset_d = multiset({'a':1, 'b':1, 'z':2})
        self.assertEqual(mset_a.symmetric_difference(mset_b, mset_c), mset_d)
        self.assertIsNot(mset_a.symmetric_difference(mset_b, mset_c), mset_a)

    def test_mul(self):
        """test multiplication"""
        mset_a = multiset({'a':1, 'z':2})
        mset_b = multiset({'a':3, 'z':6})
        self.assertEqual(mset_a*3, mset_b)
        self.assertNotEqual(mset_a*3, mset_a)

    def test_rmul(self):
        """test multiplication"""
        mset = multiset({'a':1, 'z':2})
        self.assertEqual(3*mset, mset*3)
        self.assertNotEqual(3*mset, mset)

    def test_imul(self):
        """test inplace multiplication"""
        mset_a = multiset({'a':1, 'z':2})
        mset_b = multiset({'a':3, 'z':6})
        mset_c = mset_a
        mset_a *= 3
        self.assertEqual(mset_a, mset_b)
        self.assertIs(mset_a, mset_c)

    def test_floordiv(self):
        """test division"""
        mset_a = multiset({'a':3, 'x': 2, 'z':7})
        mset_b = multiset({'a':1, 'z':2})
        self.assertEqual(mset_a//3, mset_b)
        self.assertNotEqual(mset_a//3, mset_a)

    def test_ifloordiv(self):
        """test inplace division"""
        mset_a = multiset({'a':3, 'z':7})
        mset_b = multiset({'a':1, 'z':2})
        mset_c = mset_a
        self.assertEqual(mset_a//3, mset_b)
        self.assertIs(mset_a, mset_c)


if __name__ == '__main__':
    unittest.main()
