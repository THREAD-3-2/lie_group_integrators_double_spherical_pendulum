==========================================================
 Documentation of `variable_stepsize_lie_group_integrator`
==========================================================

`variable_stepsize_lie_group_integrator <https://github.com/THREAD-3-2/variable_stepsize_lie_group_integrator>`_
is a MATLAB code for the comparison of the performance of constant and variable step size methods. 
The RKMK pair coming from Dormand–Prince method (DOPRI 5(4), also denoted as RKMK(5,4)) is compared to RKMK of order 5 (denoted by RKMK5).
The methods are tested on the N-fold 3D pendulum example. The quality of the
approximation is measured against a reference solution obtained with ODE45 from
MATLAB with a strict tolerance. The code is part of the source code developed at the `Department of Mathematical Sciences at NTNU <https://www.ntnu.edu/imf>`_, for the papers `(Celledoni, Çokaj, Leone, Murari and Owren, (2021) International Journal of Computer Mathematics) <https://doi.org/10.1080/00207160.2021.1966772>`_ and `(Celledoni, Çokaj, Leone, Murari and Owren, (2021) arXiv) <https://doi.org/10.48550/arXiv.2109.12325>`_.


Contents
========

.. toctree::
   installation
   RKMK_var_step
   tredpend
   matlab
   :maxdepth: 2


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
