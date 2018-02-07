"""Tests for the stocal.transitions module"""
import unittest
import inspect

from .abstract_test import AbstractTestCase
import stocal


class TestProcess(unittest.TestCase):
    """Process specification

    To test a custom process subclass, derive from TestProcess and
    overload the class attribute Process with the class to be tested.
    """
    Process = stocal.Process

    def test_init(self):
        """Processes can be initialized with and without reactions/rules"""
        self.Process([])
        self.Process([], [])
        self.Process([], rules=[])
        self.Process(rules=[])
        self.Process(transitions=[])

    def test_trajectory_arguments(self):
        """Process.trajectory can be called with optional arguments"""
        proc = self.Process([])
        proc.trajectory({})
        proc.trajectory({}, 1.)
        proc.trajectory({}, 1., 2.)
        proc.trajectory({}, t=1.)
        proc.trajectory({}, tmax=2.)
        proc.trajectory({}, steps=10)

    def test_trajectory_with_events(self):
        """Partly deterministic processes return an appropriate sampler"""
        proc = self.Process([stocal.Event({}, {'a':1}, 1.)])
        proc.trajectory({})


class TestRule(AbstractTestCase('Rule', stocal.Rule)):
    """Rule specification

    To test a custom rule subclass, derive from TestRule and
    overload TestRule.Rule with the rule class to be tested.
    """
    def test_init(self):
        """Supports initialization without arguments"""
        self.Rule()

    def test_infer_transitions_null(self):
        """Rule.infer_transitions infers nothing if new_species is empty."""
        if inspect.isabstract(self.Rule):
            return
        rule = self.Rule()
        self.assertFalse(list(rule.infer_transitions({}, {})))
        self.assertFalse(list(rule.infer_transitions({}, {'a': 3})))


class TestReactionRule(TestRule):
    """Test stocal.ReactionRule interface

    This abstract test case provides the test for the specification of
    ReactionRule classes. To test a concrete ReactionRule, derive a test
    case from TestReactionRule and override the class method Rule
    with the class to be tested.
    """
    Rule = stocal.ReactionRule

    def test_order_positive(self):
        """Rule.order must be greater than 0"""
        self.assertGreater(self.Rule().order, 0)

    def test_novel_reaction_inferface(self):
        """Rule.novel_reactions must not use variable argument list"""
        try:
            # python 3.5
            parameters = inspect.signature(self.Rule.novel_reactions).parameters
            self.assertFalse(any(par for par in parameters.values()
                             if par.kind == inspect.Parameter.VAR_POSITIONAL))
        except AttributeError:
            # python 2.7
            signature = inspect.getargspec(self.Rule.novel_reactions)
            self.assertIsNone(signature.varargs)

class TestTransition(AbstractTestCase('Transition', stocal.Transition)):
    """Test stocal.Transition interface

    This abstract test case provides the test for the specification of
    Transition classes. To test a concrete Transition, derive a test
    case from TestTransition and override the class method Transition
    with the class to be tested.
    """
    def test_init_positive(self):
        """Reactants and products must be positive"""
        self.Transition({'a':1}, {'z':1})
        with self.assertRaises(ValueError):
            self.Transition({'a':0}, {'z':1})
        with self.assertRaises(ValueError):
            self.Transition({'a':-1}, {'z':1})
        with self.assertRaises(ValueError):
            self.Transition({'a':1}, {'z':0})
        with self.assertRaises(ValueError):
            self.Transition({'a':1}, {'z':-1})

    def test_init_nonempty(self):
        """Either reactants or products must be specified"""
        self.Transition({'a':1}, {})
        self.Transition({'a':1}, {'z':1})
        with self.assertRaises(ValueError):
            self.Transition({}, {})

    def test_init_sequence(self):
        """Transitions can be intialized with dicts or sequences"""
        self.assertEqual(self.Transition(['a'], ['z']), self.Transition({'a':1}, {'z':1}))
        self.assertEqual(self.Transition(['a', 'a'], ['z']), self.Transition({'a':2}, {'z':1}))

    def test_equality(self):
        """Transitions are equal if their reactants and products are equal"""
        trans = self.Transition({'a':1}, {'z':1})
        self.assertEqual(trans, self.Transition({'a':1}, {'z':1}))
        self.assertNotEqual(trans, self.Transition({'b':1}, {'z':1}))
        self.assertNotEqual(trans, self.Transition({'a':1}, {'y':1}))

    def test_true_reactants(self):
        """true_reactants and true_products exclude catalysts"""
        trans1 = self.Transition({'a':1}, {'a':2})
        trans2 = self.Transition({}, {'a':1})
        self.assertNotEqual(trans1.true_reactants, trans1.reactants)
        self.assertNotEqual(trans1.true_products, trans1.products)
        self.assertEqual(trans1.true_reactants, trans2.true_reactants)
        self.assertEqual(trans1.true_products, trans2.true_products)

    def test_hash(self):
        """Equal transitions must have equal hash values"""
        trans_1 = self.Transition({'a':1}, {'z':1})
        trans_2 = self.Transition({'a':1}, {'z':1})
        self.assertEqual(hash(trans_1), hash(trans_2))


class TestEvent(TestTransition):
    """Test stocal.Event"""
    class Transition(stocal.Event):
        """Provide default values for TestTransition tests"""
        def __init__(self, reactants, products, time=1., frequency=0.):
            super(TestEvent.Transition, self).__init__(reactants, products, time, frequency)

    def test_next_occurrence_singular(self):
        """fires at given time for any earlier time and inf after"""
        time = 0.1
        event = self.Transition({}, {'a':1}, time)
        self.assertEqual(event.next_occurrence(0.), time)
        self.assertEqual(event.next_occurrence(2*time), float('inf'))

    def test_next_occurrence_periodic(self):
        """fires repeatedly with given offset and period"""
        offset = 0.1
        delta_t = 1.
        event = self.Transition({}, {'a':1}, offset, delta_t)
        for time in range(100):
            self.assertAlmostEqual(event.next_occurrence(time), time+offset)


class TestReaction(TestTransition):
    """Test stocal.Reaction interface"""
    Transition = stocal.Reaction

    def test_equality(self):
        """Reactions are equal if their propensity functions are identical"""
        class SameTransition(self.Transition):
            """Transition that should be treated equal"""
            pass
        class OtherTransition(self.Transition):
            """Transition that should be treated unequal"""
            def propensity(self, state):
                return 0.
        transition = self.Transition({'a':1}, {'z':1})
        self.assertEqual(transition, SameTransition({'a':1}, {'z':1}))
        self.assertNotEqual(transition, OtherTransition({'a':1}, {'z':1}))

    def test_propensity(self):
        """Inapplicable reactions have zero propensity"""
        transition = self.Transition({'a':1}, {'z':1})
        self.assertEqual(transition.propensity({}), 0.)

    def test_next_occurrence(self):
        """Inapplicable reactions occur in infinite time"""
        transition = self.Transition({'a':1}, {'z':1})
        self.assertEqual(transition.next_occurrence(0., {}), float('inf'))


class TestMassAction(TestReaction):
    """Test stocal.MassAction reactions"""
    class Transition(stocal.MassAction):
        """Provide default rate constant for use in TestTransition"""
        def __init__(self, reactants, products, c=1.):
            super(TestMassAction.Transition, self).__init__(reactants, products, c)

    def test_init(self):
        """rate constant must be non-negative"""
        with self.assertRaises(ValueError):
            self.Transition({'a':1}, {'z':1}, -1.)

    def test_equality_constants(self):
        """MassActions with different constants are considered different."""
        self.assertEqual(self.Transition({'a':1}, {'z':1}, 1.),
                         self.Transition({'a':1}, {'z':1}, 1.))
        self.assertNotEqual(self.Transition({'a':1}, {'z':1}, 1.),
                            self.Transition({'a':1}, {'z':1}, 2.))

    def test_equality_from_equivalent_definitions(self):
        """Assert propensity is identical for equivalent definitions

        No matter whether reactants are provided as dicts or sequences,
        the propensity of the reaction must not change.
        """
        self.assertEqual(self.Transition({'a': 2}, {}, 1.2), self.Transition(['a', 'a'], {}, 1.2))

    def test_propensity_equals(self):
        """Assert consistent propensities

        This test compares the propensities of a single reaction
        2*A --> 0 in a system with n molecules to a reaction system
        A+B --> 0, 2*A --> 0, 2*B -->0 in a system with n_a+n_b==n
        molecules. If the physical properties of the molecules are
        equal, the total propensity of the reactions has to be the same.
        """
        n_a = 12
        n_b = 8
        constant = 1.2
        a_1 = self.Transition({'a':2}, {}, constant).propensity({'a':n_a+n_b})
        a_2 = sum(r.propensity({'a':n_a, 'b':n_b}) for r in [
            self.Transition({'a':2}, {}, constant),
            self.Transition({'b':2}, {}, constant),
            self.Transition({'a':1, 'b':1}, {}, constant),
        ])
        self.assertAlmostEqual(a_1, a_2)


if __name__ == '__main__':
    unittest.main()
