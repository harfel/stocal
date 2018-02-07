# stocal tutorial
stocal is a python library for rule-based stochastic simulation (AKA [Gillespie simulations](https://en.wikipedia.org/wiki/Gillespie_algorithm)). 
Stocal works with python version 2.7 as well as version 3.5.

## Download and Installation
The latest stable release of stocal is available from the python package index:
```shell
pip install stocal
```
The development version can be obtained from github using the following commands:
```shell
git clone https://github.com/harfel/stocal.git
cd stocal
git checkout develop
sudo python setup.py install
```

## A Simple Example
We start by defining a simple stochastic process that describes the reversible dimerization of two molecules of type A into a dimer A2. Reactions are supposed to follow mass action kinetics.

To model this in stocal, we define two MassAction reactions, one that transforms two A molecules into one molecule A2 with stochastic rate constant 1.0, and one that transforms one molecule of A2 into two A molecules with stochastic rate constant 10.0. We then create a Process that groups these two reactions.
```python
from stocal import *
r1 = MassAction({'A': 2}, {'A2': 1}, 1.0)
r2 = MassAction({'A2': 1}, {'A': 2}, 10.0)
process = Process([r1, r2])
```

We can use this process, to sample stochastic trajectories. The method `Process.trajectory` instantiates a trajectory sampler for a given initial condition and stop criterion. The trajectory sampler implements the iterator protocol, so we can simply iterate through the trajectory, invoking one stochastic transition at a time. With each transition, time and state of the trajectory are properly updated: 
```python
trajectory = process.trajectory({'A':100}, steps=1000)
for transition in trajectory :
    print trajectory.time, trajectory.state['A'], trajectory.state['A2']
```
This writes out the time and state of the first 1000 steps of the stochastic trajectory as white space separated data.

## Reactions
Reactions such as MassActions above can either be invoked by specifying reactant and product stoichiometries as dictionaries, or as lists. If lists are used, the order of elements in the list is irrelevant. There is no bound to the number of molecules that can undergo a reaction. Reactants and products can also be empty to model inflow into and outflow out of the system.

Propensities of reactions are calculated by multiplying the stochastic rate constant with the number of potential reaction partners in a given state. For reactions of up to three reactants, this reads:

| Reactants | Propensity                                                  |
| --------- | ----------------------------------------------------------- |
| -         | reaction.c                                                  |
| A         | reaction.c\*state['A']                                      |
| A + B     | reaction.c\*state['A']\*state['B']                          |
| A + A     | reaction.c\*state['A']\*(state['A']-1])/2.                  |
| A + B + C | reaction.c\*state['A']\*state['B']\*state['C']              |
| A + A + B | reaction.c\*state['A']\*(state['A']-1])\*state['B']/2.      |
| A + A + A | reaction.c\*state['A']\*(state['A']-1])\*(state['A']-2])/6. |

In general, the propensity of a reaction is the stochastic rate constant times the product of the binomial coefficients to choose _n_ reaction partners, _n_ being the stoichiometry of the reactant, out of _m_ molecules, _m_ being the copy number of that reactant in the system state, for each reactant type.

## Events
stocal supports deterministic transitions that happen at a specific time, either once or periodically with a given frequency. 
```python
feed = Event([], ['A'], 0.0, 1.0)
process = Process([r1, r2, feed])
```
Here, the event `feed` will occur, and feed an A molecule to the system, at time t=0.0 and then periodically every 1.0 time units. Unlike stochastic reactions that occur with an average frequency, nondeterministic events happen at exactly the specified times.

## Rule-Based Processes
Having introduced an inflow, we next add an outflow to the model that dilutes species proportional to their abundance. We could simply add two reactions that remove molecules from the state:
```python
r4 = MassAction(['A'], [], 0.001)
r5 = MassAction(['A2'], [], 0.001)
```
However, this requires adding a dilution reaction for every chemical in model. This might become cumbersome when dealing with many species, and we might end up forgetting the dilution of one species or another.

We take this scenario as a motivation to introduce rule-based modeling. Rules can be thought of as patterns for reactions, rather than specific reactions. As such, rules generate a whole set of reactions.

Defining a rule requires to create a python [class](https://docs.python.org/2/tutorial/classes.html) with some required attributes and methods. The class needs to be derived from `ReactionRule`, which requires our subclass to have the following attributes:

| attribute | description |
| --------- | ----------- |
| `Transition` | The type of Transition that the rule generates. Here, this is the `MassAction` type |
| `novel_reactions` | A method that generates an iterable of transitions for the given reactants. |

Taking this all together, we define the following Dilution rule: 
```python
class Dilution(ReactionRule) :
    Transition = MassAction

    def novel_reactions(self, species) :
        yield self.Transition([species], [], 0.001)
```
Note the use of `yield` in the `novel_reactions` method. This [python keyword](https://docs.python.org/2/reference/simple_stmts.html#the-yield-statement) generates an on-the-fly iterable that contains all yielded items. If `yield` is unfamiliar to you, you can instead return a list of transitions without changing the behavior of the code:
```python
    def novel_reactions(self, species) :
        return [ self.Transition([species], [], 0.001) ]
```
Having defined a new rule, we can create a rule-based stochastic process by giving a second argument the Process constructor:
```python
process = Process([r1, r2, feed], [Dilution()])
```
Note here, that the second argument is a list of rule _instances_ rather than classes.

For clarity, `Process` allows its arguments to be named, and we could have written the same processes instantiation as
```python
process = Process(transitions=[r1, r2, feed], rules=[Dilution()])
```

Let us look at a more interesting case and consider a system where A molecules cannot only form dimers but polymers of any length. Any two polymers---including monomers which are really just polymers of length one---can come together to form a chain that joins these two polymers.

To model this, we define a rule class for the polymerization that generates a Polymerization reaction for any two reactants:
```python
class Polymerization(ReactionRule) :
    Transition = MassAction

    def novel_reactions(self, k, l) :
        yield self.Transition([k,l], [k+l], 10.)
```
This time, `novel_reactions` receives two reactants, `a` and `b` and yields a reaction that produces their concatenation. This way, rules can create molecular species that had not been previously in the system state!

To complete this example, we also generalize the reverse reactions and define a Hydrolysis rule that breaks a polymer at any bond. To make the model a little more interesting, we decide that the stochastic rate constants of these reactions depends on the lengths of the hydrolysis products, so that polymers are more likely to break in the middle.
```python
class Hydrolysis(ReactionRule) :
    Transition = MassAction

    def novel_reactions(self, k) :
        for i in range(1, len(k)) :
            c = 10.*i*(len(k)-i)
            yield self.Transition([k], [k[:i], k[i:]], c)
```
This time our rule employs a `for` loop to generate several reactions for each reactant---one for each potential breaking point of the polymer.

The total stochastic process, including feeding, polymerization, hybridization, and dilution is then defined by:
```python
process = Process(transitions=[feed], rules=[Dilution(), Polymerization(), Hydrolysis()])
```
Note that no change is necessary for the dilution rule, since it already generates a reaction for every chemical in the system.

## Complex States
So far, all our molecular species have been character sequences, either in the form of simple labels such as "A" and "A2", or in the form of strings. However, stocal does not require chemicals to be strings. Any immutable object can be used as a valid chemical species. Examples would be tuples, `frozensets`, or custom python classes that define a `__hash__` method and do not allow the user to alter the state of an instance. This functionality is handy when modeling chemistries that are more complex than simple molecules and polymers.

When defining custom classes to work with stocal, it is important to properly implement what is called _structural congruence_. Simply put, structurally congruent objects objects that are physically identical (congruent) even though they might differ syntactically.

As a simple example, imagine we would like to model molecular complexes, i.e. non-covalent associations of molecules. These are important, for example, in molecular biology, where many proteins form multi-protein complexes.

We could decide to model those complexes using tuples, where the tuple items correspond to the individual components of the complex. For example, the tuple `('50S', '30S')` could refer to the complexified large and small subunit that constitute the ribosome.

However, tuples are ordered sequences in python, whereas molecular complexes usually do not have a designated order of their components: `('50S', '30S')` is really just the same as `('30S', '50S')` and we need to teach this to python.

To do so, we define a custom data type (class) that provides an implementation of the equality operator `__eq__` as well as the hash function `__hash__`. Since python does not impose any semantics on custom operators, we also have to define the inequality operator `__ne__`.

The simplest way to implement structural congruence is by means of a normalization function that maps all congruent instances to an identical representation. For our molecular complexes, we could simply sort the tuple elements, thus making sure that differently ordered complexes have the same normalization.
```python
class Complex(tuple) :
    @property
    def normalized(self) :
        return tuple(sorted(self))

    def __eq__(self, other) :
        return self.normalized == other.normalized

    def __ne__(self, other) :
        return not self==other

    def __hash__(self) :
        return hash(self.normalized)
```
A complete implementation for molecular complexes might also overload the plus operator `__add__` to make sure that adding to Complexes generates a new Complex (since it would currently generate a tuple).

## Propensities
Rule-based stochastic processes bear a subtlety with regard to propensities which does not appear in regular stochastic processes. It is thus worthwhile to discuss propensity calculations in more detail.

To illustrate the issue, we extend the above polymer example to work with several types of monomers A and B, which can form polymers with mixed content, such as ABBABAA. To achieve this, we simply need to define another feed Event that provides monomers of the second type:
```python
process = Process(
    transitions=[
        Event({}, {'A': 1}, 0., 1.),
        Event({}, {'B': 1}, 0., 1.),
    ],
    rules=[Dilution(), Polymerization(), Hydrolysis()]
)
```
However, we need to decide what polymerization means and need to slightly adapt the code of our model. A (linear) polymer is a chain of interlinked polymers. Links could either be directional or non-directional. Chemical examples of directional links are ester bonds, peptide bonds, nucleic acid bonds, or any other bond where one can clearly identify a left-hand and a right-hand side in the polymer. Ether bonds, ketones and thiol bonds, on the other hand, are examples of non-directional bonds, where the molecule is rotationally symmetric along the binding site.

We have to decide whether our model features directional or non-directional polymerization. Our choice will determine which route we need to take to model polymerization accordingly.

In the case of directional bonds, two polymers _k_ and _l_ can potentially form two different polymerization products: _k+l_ and _l+k_. Therefore, the polymerization rule has to generate both reactions:
```python
class Polymerization(ReactionRule) :
    Transition = MassAction

    def novel_reactions(self, k, l) :
        yield self.Transition([k,l], [k+l], 5.)
        yield self.Transition([k,l], [l+k], 5.)
```
If _k_ and _l_ are different and if _k+l_ is different from _l+k_, this yields two reactions with propensities 5 n<sub>k</sub> n<sub>l</sub> each. If, however, _k_ equals _l_, the generated reactions are identical, each one with propensity 5/2 n<sub>k</sub> (n<sub>k</sub>-1), where the factor 1/2 comes from the binomial coefficient discussed in the section on Reactions. It is also possible for _k_ and _l_ to be different, but yet, for the reaction products _k+l_ and _l+k_ to be identical--or more precisely, structurally congruent. An example would be the molecules AB and ABAB which form the polymer ABABAB no matter which way around they bind. Since reactant and product lists in the Transition constructor are unordered lists, the two generated reactions would also be identical, each with propensity 5 n<sub>k</sub> n<sub>l</sub>.

Stocal properly detects the multiplicity of reactions, and assigns to each generated reaction a total propensity that sums up the propensities from individually generated reactions. In the example, this implies that any "left" polymer will bind any "right" polymer with about the same propensity, no matter whether the two reactants are equal or not. However, if the two possible polymerization products are indistinguishable, they will be produced with a doubled propensity.

In the case of non-directional bonds, we only have to infer the original one reaction, but we have to assert that _k+l_ and _l+k_ are structurally congruent. As we have seen before, this is best done by defining a custom type for non-directional polymers:
```python
class Polymer(str) :
    @property
    def normalized(self) :
        return min(self, ''.join(reversed(self)))
```
with the above overloads for `__eq__`, `__ne__` and `__hash__`.
The nondirectional Polymerization rule now becomes:
```python
class Polymerization(ReactionRule) :
    Transition = MassAction

    def novel_reactions(self, k, l) :
        yield self.Transition([k,l], [Polymer(k+l)], 10.)
```
In this case, propensities are calculated as in the standard Gillespie algorithm, where the propensity of a reaction with distinguishable partners is twice as big as the propensity of reactions with indistinguishable partners.

In summary, when modeling chemistries in stocal, the user does not need to bother about calculating propensities, as this is dealt with by the framework. In contrast, what the user has to pay attention to is that the textual representation of molecules properly captures the physical aspects of the modeled chemistry, i.e. define proper structural congruence relations.

## Further Documentation
The full API of stocal is available via pydoc:
```bash
pydoc stocal
```
Examples of stocal in use can be found in the stocal/examples folder.
