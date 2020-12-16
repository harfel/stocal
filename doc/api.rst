Public API
==========
.. module:: stocal

Stocal supports simulation of stochastic processes, which are either
static sets of transitions such as chemical reactions, or transitions
that are dynamically inferred from a (static) set of rules.
At the heart of stocal are thus `processes <#Processes>`_,
`transitions <#Transitions>`_, and `rules <#Rules>`_.
See the `tutotial <tutorial.html>`_ for an informal introduction to
these concepts.


Processes
---------
The Process class is the main entry point into stocal. Processes
collect transitions and rules and provide an interface to operate
over them.

.. autoclass:: Process
	:members:

State
-----

.. autoclass:: multiset
	:members:
	:special-members:

Transitions
-----------
Transitions are transformations of state elements.

.. autoclass:: Transition
	:members:

.. autoclass:: Reaction
	:members:

.. autoclass:: MassAction
	:members:

.. autoclass:: Event
	:members:


Rules
-----
.. autoclass:: Rule
	:members:

.. autoclass:: TransitionRule
	:members:


Trajectories and Ensembles
--------------------------

.. autoclass:: Trajectory
	:members:

.. py:class: Ensemble
	:members:
