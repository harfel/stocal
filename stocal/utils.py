"""Collection of utilities"""

def with_metaclass(meta, *bases):
    """Create a base class with a metaclass.

    Code taken from six (https://pypi.python.org/pypi/six).
    """
    # This requires a bit of explanation: the basic idea is to make a dummy
    # metaclass for one level of class instantiation that replaces itself with
    # the actual metaclass.
    class metaclass(meta):

        def __new__(cls, name, _, doc):
            return meta(name, bases, doc)
    return type.__new__(metaclass, 'temporary_class', (), {})
