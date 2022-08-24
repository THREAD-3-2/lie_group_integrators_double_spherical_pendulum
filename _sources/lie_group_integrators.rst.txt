
.. _lie_group_integrators:

=====================
Lie group integrators
=====================

.. _diff_eqs_in_manifolds:

The formulation of differential equations on manifolds
------------------------------------------------------


Lie group integrators solve differential equations whose solution evolve on a manifold :math:`\mathcal{M}`. 
This means that we seek a curve :math:`y(t)\in\mathcal{M}` whose tangent at any point coincides with a vector field :math:`F\in\mathcal{X}(\mathcal{M})` and passing through a designated initial value
:math:`y_0` at :math:`t=t_0`

.. math::
    :name: eq:1

    \begin{align}
        \dot{y}(t) = F|_{y(t)},\qquad y(t_0)=y_0.    
    \end{align}

Before addressing numerical methods for solving :ref:`(1) <eq:1>` it is necessary to introduce a convenient way of representing the vector field :math:`F`. 
We furnish :math:`\mathcal{M}` with a transitive action :math:`\Psi: G \times \mathcal{M} \rightarrow \mathcal{M}` by some Lie group :math:`G` of dimension :math:`d\geq\dim \mathcal{M}`. 
We denote the action of :math:`g` on :math:`m` as :math:`g\cdot m`, i.e. :math:`g\cdot m=\Psi(g,m)`.
Let :math:`\mathfrak{g}` be the Lie algebra of :math:`G`, and denote by :math:`\exp: \mathfrak{g}\rightarrow G` the exponential map. 
We define  :math:`\Psi_*:\mathfrak{g}\rightarrow\mathcal{X}(\mathcal{M})` to be the infinitesimal generator of the action, i.e.

.. math::
    :name: eq:

    \begin{align}
        \left.F_\xi\right|_m=  \left.\Psi_*(\xi)\right|_m = \left.\frac{d}{dt}\right|_{t=0} \Psi(\exp(t\xi), m)
    \end{align}

The transitivity of the action now ensures that :math:`\left.\Psi_*(\mathfrak{g})\right|_m=T_m\mathcal{M}` for any :math:`m\in\mathcal{M}`, such that any tangent vector :math:`v_m\in T_m\mathcal{M}` can be represented as :math:`v_m=\left.\Psi_*(\xi_v)\right|_m` for some :math:`\xi_v\in\mathfrak{g}` (:math:`\xi_v` may not be unique). 
Consequently, for any vector field :math:`F\in\mathcal{X}(\mathcal{M})` there exists a map :math:`f:\mathcal{M}\rightarrow\mathfrak{g}`
}} such that

.. math::
    :name: eq:

    \begin{align}
        F|_m = \left.\Psi_*(f(m))\right|_m,\quad\text{for all}\; m\in \mathcal{M}
    \end{align}

This is the original tool from `(H. Munthe-Kaas, (1999) Appl. Numer. Math 29) <https://doi.org/10.1016/S0168-9274(98)00030-0>`_ for representing a vector field on a manifold with a group action.