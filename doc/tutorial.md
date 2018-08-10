# stocal tutorial

stocal is a python library for rule-based stochastic simulation (AKA
[Gillespie simulations](https://en.wikipedia.org/wiki/Gillespie_algorithm)). 
Stocal works with python version 2.7 as well as version 3.5.

## Download and Installation

The latest stable release of stocal is available from the python package
index:

```shell
pip install stocal
```

The development version can be obtained from github using the following
commands:

```shell
git clone https://github.com/harfel/stocal.git
cd stocal
git checkout develop
sudo python setup.py install
```

## A Simple Example

We start by defining a simple stochastic process that describes the
reversible dimerization of two molecules A into a dimer A2. Reactions
are supposed to follow mass action kinetics.

To model this in stocal, we define two MassAction reactions, one that
transforms two A molecules into one molecule A2 with stochastic rate
constant 1.0, and one that transforms one molecule of A2 into two A
molecules with stochastic rate constant 10.0. We then create a Process
that groups these two reactions.

```python
from stocal import *
r1 = MassAction({'A': 2}, {'A2': 1}, 1.0)
r2 = MassAction({'A2': 1}, {'A': 2}, 10.0)
process = Process([r1, r2])
```

We can use this process, to sample stochastic trajectories. The method
`Process.sample` instantiates a trajectory sampler for a given
initial condition and stop criterion. The trajectory sampler implements
the iterator protocol, so we can simply iterate through the trajectory,
invoking one stochastic transition at a time. With each transition,
time and state of the trajectory are properly updated: 

```python
trajectory = process.sample({'A':100}, steps=1000)
for dt, transitions in trajectory :
    print trajectory.time, trajectory.state['A'], trajectory.state['A2']
```

This writes out the time and state of the first 1000 steps of the
stochastic trajectory as white space separated data.

## Reactions

Reactions such as MassActions above can either be invoked by specifying
reactant and product stoichiometries as dictionaries, or as lists.
If lists are used, the order of elements in the list is irrelevant.
There is no bound to the number of molecules that can undergo a
reaction. Reactants and products can also be empty to model inflow into
and outflow out of the system.

Propensities of reactions are calculated by multiplying the stochastic
rate constant with the number of potential reaction partners in a given
state. For reactions of up to three reactants, this reads:

| Reactants | Propensity                                                  |
| --------- | ----------------------------------------------------------- |
| -         | reaction.c                                                  |
| A         | reaction.c\*state['A']                                      |
| A + B     | reaction.c\*state['A']\*state['B']                          |
| A + A     | reaction.c\*state['A']\*(state['A']-1])/2.                  |
| A + B + C | reaction.c\*state['A']\*state['B']\*state['C']              |
| A + A + B | reaction.c\*state['A']\*(state['A']-1])\*state['B']/2.      |
| A + A + A | reaction.c\*state['A']\*(state['A']-1])\*(state['A']-2])/6. |

In general, the propensity of a reaction is the stochastic rate constant
times the product of the binomial coefficients to choose _n_ reaction
partners, _n_ being the stoichiometry of the reactant, out of _m_
molecules, _m_ being the copy number of that reactant in the system
state, for each reactant species.

## Events

Stocal supports deterministic transitions that happen at a specific
time, either once or periodically with a given frequency. 

```python
feed = Event([], ['A'], 0.0, 1.0)
process = Process([r1, r2, feed])
```

Here, the event `feed` will occur, and feed an A molecule to the system,
at time t=0.0 and then periodically every 1.0 time units. Unlike
stochastic reactions that occur with an average frequency,
nondeterministic events happen at exactly the specified times.

## Rule-Based Processes

Having introduced an inflow, we next add an outflow to the model that
dilutes species proportional to their abundance. We could simply add
two reactions that remove molecules from the state:

```python
r4 = MassAction(['A'], [], 0.001)
r5 = MassAction(['A2'], [], 0.001)
```

However, this requires adding a dilution reaction for every chemical in
model. This might become cumbersome when dealing with many species, and
we might end up forgetting the dilution of one species or another.

We take this scenario as a motivation to introduce rule-based modeling.
Rules can be thought of as patterns for reactions, rather than specific
reactions. As such, rules generate a whole set of reactions.

Defining a rule requires to create a python
[class](https://docs.python.org/2/tutorial/classes.html) with some
required attributes and methods. The class needs to be derived from
`TransitionRule`, which requires our subclass to have the following
attributes:

| attribute         | description                                                                         |
| ----------------- | ----------------------------------------------------------------------------------- |
| `Transition`      | The type of Transition that the rule generates. Here, this is the `MassAction` type |
| `novel_reactions` | A method that generates an iterable of transitions for the given reactants.         |

Taking this all together, we define the following Dilution rule: 

```python
class Dilution(TransitionRule) :
    Transition = MassAction

    def novel_reactions(self, species) :
        yield self.Transition([species], [], 0.001)
```

Note the use of `yield` in the `novel_reactions` method. This
[python keyword](https://docs.python.org/2/reference/simple_stmts.html#the-yield-statement)
generates an on-the-fly iterable that contains all yielded items. If
`yield` is unfamiliar to you, you can instead return a list of
transitions without changing the behavior of the code:

```python
    def novel_reactions(self, species) :
        return [ self.Transition([species], [], 0.001) ]
```

*New in version 1.1:* In python3, the transition type of a rule can
alternatively be provided as return type annotation of the
`novel_reactions` method. For example:

```python
from typing import Iterator

class Dilution(TransitionRule) :
    def novel_reactions(self, species) -> Iterator[MassAction]:
        yield MassAction([species], [], 0.001)
```

Having defined a new rule, we can create a rule-based stochastic process
by giving a second argument to the Process constructor:

```python
process = Process([r1, r2, feed], [Dilution()])
```

Note here, that the second argument is a list of rule _instances_ rather
than classes.

For clarity, `Process` allows its arguments to be named, and we could
have written the same process instantiation as

```python
process = Process(transitions=[r1, r2, feed], rules=[Dilution()])
```

Let us look at a more interesting case and consider a system where A
molecules cannot only form dimers but polymers of any length. Any two
polymers---including monomers which are really just polymers of length
one---can come together to form a chain that joins these two polymers.

To model this, we define a rule class for the polymerization that
generates a Polymerization reaction for any two reactants:

```python
class Polymerization(TransitionRule) :
    Transition = MassAction

    def novel_reactions(self, k, l) :
        yield self.Transition([k,l], [k+l], 10.)
```

This time, `novel_reactions` receives two reactants, `k` and `l` and
yields a reaction that produces their concatenation. This way, rules
can create molecular species that had not been previously in the system
state!

To complete this example, we also generalize the reverse reaction and
define a Hydrolysis rule that breaks a polymer at any bond. To make the
model a little more interesting, we decide that the stochastic rate
constants of these reactions depends on the lengths of the hydrolysis
products, so that polymers are more likely to break in the middle.

```python
class Hydrolysis(TransitionRule) :
    Transition = MassAction

    def novel_reactions(self, k) :
        for i in range(1, len(k)) :
            c = 10.*i*(len(k)-i)
            yield self.Transition([k], [k[:i], k[i:]], c)
```

This time our rule employs a `for` loop to generate several reactions
for each reactant---one for each potential breaking point of the
polymer.

The total stochastic process, including feeding, polymerization,
hybridization, and dilution is then defined by:

```python
process = Process(transitions=[feed],
                  rules=[Dilution(), Polymerization(), Hydrolysis()])
```

Note that no change is necessary for the dilution rule, since it already
generates a reaction for every chemical in the system.

*New in version 1.2:* Rule-based processes that expand into a finite
set of transitions can be flattened into equivalent static processes
that employ specific transitions rather than general rules:

```python
process = Process(rules=[Dilution()])
flat_process = process.flatten(['a', 'b', 'c'])
```
This will generate a new process objects where the original rule is
expanded into three transitions, each one modelling the specific
dilution of one of the provided molecular species.


## Complex States

So far, all our molecular species have been character sequences, either
in the form of simple labels such as "A" and "A2", or in the form of
strings. However, stocal does not require chemicals to be strings. Any
immutable object can be used as a valid chemical species. Examples would
be tuples, `frozensets`, or custom python classes that define a
`__hash__` method and do not allow the user to alter the state of an
instance. This functionality is handy when modeling chemistries that are
more complex than simple molecules and polymers.

When defining custom classes to work with stocal, it is important to
properly implement what is called _structural congruence_. Simply put,
structurally congruent objects objects that are physically identical
(congruent) even though they might differ syntactically.

As a simple example, imagine we would like to model molecular complexes,
i.e. non-covalent associations of molecules. These are important, for
example, in molecular biology, where many proteins form multi-protein
complexes.

We could decide to model those complexes using tuples, where the tuple
items correspond to the individual components of the complex. For
example, the tuple `('50S', '30S')` could refer to the complexified
large and small subunit that constitute the ribosome.

However, tuples are ordered sequences in python, whereas molecular
complexes usually do not have a designated order of their components:
`('50S', '30S')` is really just the same as `('30S', '50S')` and we need
to teach this to python.

To do so, we define a custom data type (class) that provides an
implementation of the equality operator `__eq__` as well as the hash
function `__hash__`. Since python does not impose any semantics on
custom operators, we also have to define the inequality operator
`__ne__`.

The simplest way to implement structural congruence is by means of a
normalization function that maps all congruent instances to an identical
representation. For our molecular complexes, we could simply sort the
tuple elements, thus making sure that differently ordered complexes have
the same normalization.

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

A complete implementation for molecular complexes might also overload
he plus operator `__add__` to make sure that adding to Complexes
generates a new Complex (since it would currently generate a tuple).

## Propensities

Rule-based stochastic processes bear a subtlety with regard to
propensities which does not appear in regular stochastic processes. It
s thus worthwhile to discuss propensity calculations in more detail.

To illustrate the issue, we extend the above polymer example to work
with several types of monomers A and B, which can form polymers with
mixed content, such as ABBABAA. To achieve this, we simply need to
define another feed Event that provides monomers of the second type:

```python
process = Process(
    transitions=[
        Event({}, {'A': 1}, 0., 1.),
        Event({}, {'B': 1}, 0., 1.),
    ],
    rules=[Dilution(), Polymerization(), Hydrolysis()]
)
```

However, we need to decide what polymerization means and need to
slightly adapt the code of our model. A (linear) polymer is a chain of
interlinked monomers. Links could either be directional or
non-directional. Chemical examples of directional links are ester bonds,
peptide bonds, nucleic acid bonds, or any other bond where one can
clearly identify a left-hand and a right-hand side in the polymer. Ether
bonds, ketones and thiol bonds, on the other hand, are examples of
non-directional bonds, where the molecule is rotationally symmetric
along the binding site.

We have to decide whether our model features directional or
non-directional polymerization. Our choice will determine which route we
need to take to model polymerization appropriately.

In the case of directional bonds, two polymers _k_ and _l_ can
potentially form two different polymerization products: _k+l_ and _l+k_.
Therefore, the polymerization rule has to generate both reactions:

```python
class Polymerization(TransitionRule) :
    Transition = MassAction

    def novel_reactions(self, k, l) :
        yield self.Transition([k,l], [k+l], 5.)
        yield self.Transition([k,l], [l+k], 5.)
```

If _k_ and _l_ are different and if _k+l_ is different from _l+k_, this
yields two reactions with propensities 5 n<sub>k</sub> n<sub>l</sub>
each. If, however, _k_ equals _l_, the generated reactions are
identical, each one with propensity 5/2 n<sub>k</sub> (n<sub>k</sub>-1),
where the factor 1/2 comes from the binomial coefficient discussed in
the section on Reactions. It is also possible for _k_ and _l_ to be
different, but yet, for the reaction products _k+l_ and _l+k_ to be
identical--or more precisely, structurally congruent. An example would
be the molecules AB and ABAB which form the polymer ABABAB no matter
which way around they bind. Since reactant and product lists in the
Transition constructor are unordered lists, the two generated reactions
would also be identical, each with propensity
5 n<sub>k</sub> n<sub>l</sub>.

Stocal properly detects the multiplicity of reactions, and assigns to
each generated reaction a total propensity that sums up the propensities
from individually generated reactions. In the example, this implies that
any "left" polymer will bind any "right" polymer with about the same
propensity, no matter whether the two reactants are equal or not.
However, if the two possible polymerization products are
indistinguishable, they will be produced with a doubled propensity.

In the case of non-directional bonds, we only have to infer the original
one reaction, but we have to assert that _k+l_ and _l+k_ are
structurally congruent. As we have seen before, this is best done by
defining a custom type for non-directional polymers:

```python
class Polymer(str) :
    @property
    def normalized(self) :
        return min(self, ''.join(reversed(self)))
```

with the above overloads for `__eq__`, `__ne__` and `__hash__`.
The nondirectional Polymerization rule now becomes:

```python
class Polymerization(TransitionRule) :
    Transition = MassAction

    def novel_reactions(self, k, l) :
        yield self.Transition([k,l], [Polymer(k+l)], 10.)
```

In this case, propensities are calculated as in the standard Gillespie
algorithm, where the propensity of a reaction with distinguishable
partners is twice as big as the propensity of reactions with
indistinguishable partners.

In summary, when modeling chemistries in stocal, the user does not need
to bother about calculating propensities, as this is dealt with by the
framework. In contrast, what the user has to pay attention to is that
the textual representation of molecules properly captures the physical
aspects of the modeled chemistry, i.e. define proper structural
congruence relations.

## Typed Reactions

*New in version 1.1.*

If you want reaction rules to only generate reactions among certain
types of molecular species, stocal supports molecular types and typed
reaction rules. For this example, we look into modelling the association
of proteins with mRNA's. We want to define a rule for the association of
an arbitrary protein with an arbitrary mRNA. 

With the above TransitionRule's we would need to constantly check whether
the species supplied to `TransitionRule.novel_reactions` are indeed
proteins and RNA's and only yield a transition in case they are. Not
knowing which argument of the reactant combination is the protein and
which the RNA further complicates the code.

```python
class Association(TransitionRule):
    Transition = MassAction

    def novel_reactions(self, k, l):
        if is_protein(k) and is_rna(l):
            yield self.Transition([k, l], (k,l), 1.)
        elif is_rna(k) and is_protein(l):
            yield self.Transition([k, l], (l,k), 1.)
```

For these common situations, stocal offers species types and typed
rules. In stocal, the type of a species is simply its python type. So
far, we have encountered species typed as strings, Complexes, and
Polymers. Here, we define two molecule types `Protein` and `Rna` which
are simply subclasses of `str`:

```python
class Protein(str):
    pass

class Rna(str):
    pass
```

We can now write a typed `TransitionRule` for their association, simply
by setting the optional TransitionRule attribute `signature` to the list
of types that the rule should accept. When defining a signature, it must
have the same number of elements as the `novel_reactions` method.
`novel_reactions` will now only be called with arguments that adhere to
the type given in the signature. In our case, writing the rule becomes
as simple as:

```python
class Association(TransitionRule):
    Transition = MassAction
    signature = [Protein, Rna]

    def novel_reactions(self, protein, rna):
        yield self.Transition([protein, rna], [(protein,rna)], 1.)
```

In python3, type annotations can alternatively be used to specify the
rule signature:

```python
from typing import Iterator

class Association(TransitionRule):
    def novel_reactions(self, protein: Protein, rna: Rna) -> Iterator[MassAction]:
        yield MassAction([protein, rna], [(protein,rna)], 1.)
```

## Reactions with time-dependent reaction rates

*New in version 1.1.*

Stocal supports the definition of reactions with time-dependent reaction
rates. Reactions of this kind appear naturally when  parameters of the
reaction environment such as temperature or volume change over time.

As an example, let us consider reactions taking placing in a linearly
expanding reaction vessel:
```python
def volume(time, V0=1.0, dV=0.1):
    return V0 + dV*time
```

We can now define a volume dependent variant of mass action reactions.
To do so, we define a subclass of `MassAction` that overloads the
`propensity` method.
For autonomous reactions (those whose reaction rate constant is
independent of time), this function takes the `state` as sole argument.
For non-autonomous reactions, the signature is expanded to take the
`time` as a second argument. The `propensity` method can then make
use of the time argument however the user sees fit.

To accurately model volume dependency, we need to divide the reaction
rate by the volume for all but one of the reactions, i.e.
unimolecular reactions are volume independent, bimolecular reactions
are inversely proportional to the volume, a.s.o.

Taking it all together, the Reaction class reads:
```python
class VolumeDependentMassAction(MassAction):
    def propensity(self, state, time):
        a = super(VolumeDependentMassAction, self).propensity(state)
        order = sum(self.reactants.values())
        return a / volume(time)**(order-1)
```

We can now use `VolumeDependentMassAction` in any place where we
have used default `MassAction` reactions before.
stocal/examples/temperature_cycle.py gives an example of how reactions
can be modified to take changing temperature instead of volumes instead.

## Stochastic simulation algorithms

stocal ships with several variants of the stochastic simulation algorithm,
refered to as sampler. A call to `Process.sample` inspects the
underlying process and will instantiate an appropriate sampler.
Currently, this creates an instance of Gibson and Bruck's next reaction
method, unless at least one transition of the process is time-dependent
(in which case the method creates an instance of Anderon's method).

If you want to control which simulation algorithm is instantiated, you
can instantiate the desired sampler directly, as in, e.g.,

```python
sampler = algorithms.DirectMethod(process, state, tmax=100.)
for dt, transitions in sampler:
    print(dt, transitions)
```

Currently, stocal provides the following samplers:

| algorithm          | description                                                                                                    |
|------------------- | -------------------------------------------------------------------------------------------------------------- |
|DirectMethod        | Original Gillespie algorithm                                                                                   |
|FirstReactionMethod | Stochastic simulation algorithm that can operate account for scheduled events                                  |
|NextReactionMethod  | Variant of FirstReactionMethod with improved performance *(new in version 1.2)*                                |
|AndersonMethod      | Variant of NextReactionMethod that allows for propensity functions to be time-dependent *(new in version 1.1)* |
|CaoMethod           | An (inexact) tau-leaping variant of SSA -- available in stocal.experimental.tauleap *(new in version 1.2)*     |

Please refer to the class documentation for information about the exact
implementation and reference publication.

If you want to implement your own stochastic simulation algorithm, it
should be programmed against the interface defined by
`stocal.algorithms.StochasticSimulationAlgorithm`.

## Further Documentation

The full API of stocal is available via pydoc:

```bash
pydoc stocal
```

Examples of stocal in use can be found in the stocal/examples folder.
