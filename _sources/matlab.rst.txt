.. _matlab:

=============
 MATLAB code
=============

This MATLAB code is documented with Sphinx
using the `matlabdomain extension <https://github.com/sphinx-contrib/matlabdomain/blob/master/README.rst>`_.
 
 
src
========

|

.. mat:automodule:: src

:mod:`src` module contains the following source code files:
    
.. mat:autoscript:: src.main

|


src/integrators
===============

|

.. mat:autofunction:: src.integrators.RK4

|

.. mat:autofunction:: src.integrators.LieEuler

|

.. mat:autofunction:: src.integrators.LieEulerHeun

|

.. mat:autofunction:: src.integrators.CommFreeRKMK4 

|

.. mat:autofunction:: src.integrators.RKMK4 

|

.. mat:autofunction:: src.integrators.TwoCommRKMK4 

|

.. mat:autofunction:: src.integrators.RKMK3 

|

src/Lie_group_functions
=======================

|

.. mat:autofunction:: src.Lie_group_functions.expRodrigues

|

.. mat:autofunction:: src.Lie_group_functions.expSE3

|

.. mat:autofunction:: src.Lie_group_functions.exponentialSE3N

|

.. mat:autofunction:: src.Lie_group_functions.actionSE3

|

.. mat:autofunction:: src.Lie_group_functions.actionSE3N

|

.. mat:autofunction:: src.Lie_group_functions.dexpinvSE3

|

.. mat:autofunction:: src.Lie_group_functions.dexpinvSE3N

|

.. mat:autofunction:: src.Lie_group_functions.commutatorSE3

|

.. mat:autofunction:: src.Lie_group_functions.commutatorSE3N

|

src/equations_of_motion
=======================

|

.. mat:autofunction:: src.equations_of_motion.fManiToAlgebra

|

.. mat:autofunction:: src.equations_of_motion.assembleF

|

.. mat:autofunction:: src.equations_of_motion.assembleM

|

.. mat:autofunction:: src.equations_of_motion.assembleR

|

.. mat:autofunction:: src.equations_of_motion.FuncQ

|

.. mat:autofunction:: src.equations_of_motion.FuncW

|

.. mat:autofunction:: src.equations_of_motion.initializeStat

|

.. mat:autofunction:: src.equations_of_motion.initializeSE3N

|

.. mat:autofunction:: src.equations_of_motion.potential

|

src/experiments
===============

|

.. mat:autofunction:: src.experiments.checkConvergenceRate

|

.. mat:autofunction:: src.experiments.checkTangency

|

.. mat:autofunction:: src.experiments.tangentBehaviour

|

.. mat:autofunction:: src.experiments.compareNorms

|

src/helpful_functions
=====================

|

.. mat:autofunction:: src.helpful_functions.extractq

|

.. mat:autofunction:: src.helpful_functions.extractw

|

.. mat:autofunction:: src.helpful_functions.hat

|

.. mat:autofunction:: src.helpful_functions.getNorms

|

.. mat:autofunction:: src.helpful_functions.getVec

|

.. mat:autofunction:: src.helpful_functions.getBlock

|

.. mat:autofunction:: src.helpful_functions.reorder
