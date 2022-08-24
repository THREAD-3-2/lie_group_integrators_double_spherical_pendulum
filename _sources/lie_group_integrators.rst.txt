
.. _lie_group_integrators:

=====================
Lie group integrators
=====================

.. _diff_eqs_in_manifolds:

The formulation of differential equations on manifolds
------------------------------------------------------


Lie group integrators solve differential equations whose solution evolve on a manifold :math:`\mathcal{M}`. For ease of notation we restrict the discussion to the case of autonomous vector fields, although allowing for explicit :math:`t`-dependence could easily have been included.
This means that we seek a curve :math:`y(t)\in\mathcal{M}` whose tangent at any point coincides with a vector field :math:`F\in\mathcal{X}(\mathcal{M})` and passing through a designated initial value
:math:`y_0` at :math:`t=t_0`

.. math::
    :name: eq:11

    \begin{align}
        \dot{y}(t) = F|_{y(t)},\qquad y(t_0)=y_0.    
    \end{align}

Before addressing numerical methods for solving :ref:`<eq:11>` it is necessary to introduce a convenient way of representing the vector field :math:`F`. There are different ways of doing this. 
%One is to furnish :math:`\mathcal{M}` with a transitive action by some Lie group :math:`G` of dimension :math:`d\geq\dim \mathcal{M}`.
One is to furnish :math:`\mathcal{M}` with a transitive action :math:`\Psi: G \times \mathcal{M} \rightarrow \mathcal{M}` by some Lie group :math:`G` of dimension :math:`d\geq\dim \mathcal{M}`. We denote the action of :math:`g` on :math:`m` as :math:`g\cdot m`, i.e. :math:`g\cdot m=\Psi(g,m)`.
Let :math:`\mathfrak{g}` be the Lie algebra of :math:`G`, and denote by :math:`\exp: \mathfrak{g}\rightarrow G` the exponential map. We define  :math:`\Psi_*:\mathfrak{g}\rightarrow\mathcal{X}(\mathcal{M})` to be the infinitesimal generator of the action, i.e.

.. math::
    :name: eq:

    \begin{align}
        \left.F_\xi\right|_m=  \left.\Psi_*(\xi)\right|_m = \left.\frac{d}{dt}\right|_{t=0} \Psi(\exp(t\xi), m)
    \end{align}

The transitivity of the action now ensures that :math:`\left.\Psi_*(\mathfrak{g})\right|_m=T_m\mathcal{M}` for any :math:`m\in\mathcal{M}`, such that any tangent vector :math:`v_m\in T_m\mathcal{M}` can be represented as :math:`v_m=\left.\Psi_*(\xi_v)\right|_m` for some :math:`\xi_v\in\mathfrak{g}` (:math:`\xi_v` may not be unique). Consequently, for any vector field :math:`F\in\mathcal{X}(\mathcal{M})` there exists a map :math:`f:\mathcal{M}\rightarrow\mathfrak{g}`
}} such that

.. math::
    :name: eq:

    \begin{align}
        F|_m = \left.\Psi_*(f(m))\right|_m,\quad\text{for all}\; m\in \mathcal{M}
    \end{align}


This is the original tool \cite{munthe-kaas99hor} for representing a vector field on a manifold with a group action.
Another approach was used in \cite{crouch93nio} where a set of {\em frame vector fields} :math:`E_1,\ldots, E_d` in :math:`\mathcal{X}(\mathcal{M})` was introduced 
assuming that for every :math:`m\in \mathcal{M}`, 

.. math::
    :name: eq:

    \begin{align}
        \text{span}\{\left.E_1\right|_m,\ldots,\left.E_d\right|_m\}= T_m \mathcal{M}.
    \end{align}

Then, for any vector field :math:`F\in\mathcal{X}(\mathcal{M})` there are, in general non-unique, functions :math:`f_i:\mathcal{M}\rightarrow \mathbb{R}`, 
\RE{which can be chosen with the same regularity as :math:`F`,} such that

.. math::
    :name: eq:

    \begin{align}
        F|_m = \sum_{i=1}^d f_i(m) \left.E_i\right|_m.
    \end{align}

A fixed vector :math:`\xi\in\mathbb{R}^d` will define a vector field :math:`F_\xi` on :math:`\mathcal{M}` similar to \eqref{frozenaction}

.. math::
    :name: eq:1

    \begin{align}
        \left.F_{\xi}\right|_m = \sum_{i=1}^d \xi_i E_i|_m
    \end{align}

If :math:`\xi_i=f_i(p)` for some :math:`p\in\mathcal{M}`, the corresponding :math:`F_\xi` will be a vector field in the linear span of the frame which coincides with :math:`F` at the point :math:`p`. Such a vector field was named by \cite{crouch93nio} as a \emph{the vector field frozen at :math:`p`}.

The two formulations just presented are in many cases connected, and can then be used in an equivalent manner.
Suppose that :math:`e_1,\ldots,e_d` is a basis of the Lie algebra :math:`\mathfrak{g}`, then we can simply define frame vector fields as
:math:`E_i = \Psi_*(e_i)` and the vector field we aim to describe is, 

.. math::
    :name: eq:1

    \begin{align}
        F|_m=\left.\Psi_*(f(m))\right|_m= \left.\Psi_*(\sum_i f_i(m)e_i)\right|_m=\sum_i f_i \left.E_i\right|_m.
    \end{align}

As mentioned above there is a non-uniqueness issue when defining a vector field by means of a group action or a frame.
A more fundamental description can be obtained using the machinery of connections. The assumption is that the simply connected manifold :math:`\mathcal{M}` is equipped with a connection which is flat and has constant torsion. 
Then :math:`F_p`, the frozen vector field of :math:`F` at :math:`p` defined above, can be defined as the unique element :math:`F_p\in\mathcal{X}(\mathcal{M})` satisfying
\begin{enumerate}
\item :math:`F_p|_p=F|_p`
\item :math:`\nabla_X F_p=0` for any :math:`X\in\mathcal{X}(\mathcal{M})`.
\end{enumerate}
So :math:`F_p` is the vector field that coincides with :math:`F` at :math:`p` and is parallel transported to any other point on :math:`\mathcal{M}` by the connection :math:`\nabla`. Since the connection is flat, the parallel transport from the point :math:`p` to another point :math:`m\in\mathcal{M}` does not depend on the chosen path between the two points.
For further details, see e.g.
\cite{lundervold15oas}.