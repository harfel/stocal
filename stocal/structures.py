"""Collection of useful data structures"""

from collections import Mapping
from numbers import Number

class multiset(dict):
    """A multiset implementation

    Provides a feature rich implementation of multisets.
    (c.f. https://en.wikipedia.org/wiki/Multiset)
    """
    def __init__(self, *args, **opts):
        if len(args) > 1:
            raise TypeError("multiset expects at most 1 argument, got %d" % len(args))
        arg = args[0] if args else {}
        if isinstance(arg, Mapping):
            super(multiset, self).__init__(arg, **opts)
        else:
            items = {item:arg.count(item) for item in set(arg)}
            super(multiset, self).__init__(items, **opts)
        if not all(isinstance(value, Number) for value in self.values()):
            raise TypeError("multiset values must be numbers.")

    @property
    def domain(self):
        """The underlying domain (set)"""
        return set(self)

    def __repr__(self):
        return '%s(%s)' % (type(self).__name__, super(multiset, self).__repr__())

    def __eq__(self, other):
        return (all(other[item] == count for item, count in self.items())
                and all(self[item] == count for item, count in other.items()))

    def __ne__(self, other):
        return not self.__eq__(other)

    def __le__(self, other):
        return all(other[item] >= count for item, count in self.items())

    def __lt__(self, other):
        return self <= other and self != other

    def __ge__(self, other):
        return other <= self

    def __gt__(self, other):
        return other < self

    def __len__(self):
        return sum(self.values())

    def __contains__(self, item):
        if isinstance(item, multiset):
            return item <= self
        else:
            return super(multiset, self).__contains__(item)

    def __getitem__(self, item):
        return self.get(item, 0)

    def __setitem__(self, item, count):
        if count:
            super(multiset, self).__setitem__(item, count)
        else:
            del self[item]

    def __delitem__(self, item):
        if item in self:
            super(multiset, self).__delitem__(item)

    def __add__(self, other):
        result = type(self)(self)
        result += other
        return result

    def __iadd__(self, other):
        for item, count in other.items():
            self[item] += count
        return self

    def __sub__(self, other):
        result = type(self)(self)
        result -= other
        return result

    def __isub__(self, other):
        for item, count in other.items():
            if self[item] > count:
                self[item] -= count
            else:
                del self[item]
        return self

    def __mul__(self, factor):
        result = type(self)(self)
        result *= factor
        return result

    def __imul__(self, factor):
        if not factor:
            self.clear()
        else:
            for item in self:
                self[item] *= factor
        return self

    def __rmul__(self, factor):
        result = type(self)(self)
        result *= factor
        return result

    def __floordiv__(self, factor):
        result = type(self)(self)
        result //= factor
        return result

    def __ifloordiv__(self, factor):
        for item in set(self):
            self[item] //= factor
            if not self[item]:
                del self[item]
        return self

    def update(self, *args, **opts):
        """update multiset with values from mappings"""
        super(multiset, self).update(*args, **opts)
        deletes = [item for item, count in self.items() if not count]
        for item in deletes:
            del self[item]

    def union(self, *mappings):
        """The union (sum) of self and all other multisets"""
        result = type(self)(self)
        for mapping in mappings:
            result += mapping
        return result

    def difference(self, *mappings):
        """The difference of self and all other multisets"""
        result = type(self)(self)
        for mapping in mappings:
            result -= mapping
        return result

    def symmetric_difference(self, *mappings):
        """The symmetric difference of self and all other multisets"""
        result = type(self)(self)
        for mapping in mappings:
            for item, count in mapping.items():
                if result[item] > count:
                    result[item] -= count
                elif result[item] < count:
                    result[item] = count - result[item]
                else:
                    del result[item]
        return result
