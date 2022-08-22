.. _tredpend:

======================
The N-fold 3D pendulum
======================

We describe here the specific problem of a chain of :math:`N` connected 3D pendulums, whose dynamics evolves on :math:`(TS^2)^N`.
The dynamics of this mechanical system is described in terms of a Lie group :math:`G` acting transitively on the phase space :math:`\mathcal{M}`. 
The equations of motion are presented in terms of the infinitesimal generator of the transitive action.

.. _eom:

Equations of motion
-------------------

Let us consider a  chain of :math:`N` pendulums subject to constant gravity :math:`g`. The system is modeled by :math:`N` rigid, massless links serially connected by spherical joints, with the first link connected to a fixed point placed at the origin of the ambient space :math:`\mathbb{R}^3`. We neglect friction and interactions among the pendulums. 

The modeling part comes from `(Lee, Leok and McClamroch, (2018)) <https://doi.org/10.1007/978-3-319-56953-6>`_ and we omit details. We denote by :math:`q_i\in S^2` the configuration vector of the :math:`i-th` mass, :math:`m_i`, of the chain. Following `(Lee, Leok and McClamroch, (2018)) <https://doi.org/10.1007/978-3-319-56953-6>`_, we express the Euler–Lagrange equations for our system in terms of the configuration variables :math:`(q_1,\dots,q_N)\in (S^2)^N\subset\mathbb{R}^{3N}`, and their angular velocities :math:`(\omega_1,...,\omega_N)\in T_{q_1}S^2\times ... \times T_{q_N}S^2\subset\mathbb{R}^{3N}`, defined be the following kinematic equations:

.. math::
    :name: eq:1 
    
    \begin{align}
        \dot{q}_i = \omega_i\times q_i, \quad i=1,\dots,N.
    \end{align}

The Euler–Lagrange equations of the system can be written as

.. math::
    :name: eq:2
    
    \begin{align}
        R(q)\dot{\omega} = \left[\sum_{\substack{j=1\\ j\neq i}}^N M_{ij}|\omega_j|^2\hat{q}_i q_j - \Big(\sum_{j=i}^N m_j\Big)gL_i \hat{q}_i e_3 \right]_{i=1,...,N} = \begin{bmatrix}r_1\\ \vdots \\ r_N \end{bmatrix}\in\mathbb{R}^{3N},
    \end{align}
    
where :math:`R(q)\in\mathbb{R}^{3N\times 3N}` is a symmetric block matrix defined as

.. math::

    \begin{align}
        R(q)_{ii} = \Big(\sum_{j=i}^Nm_j\Big)L_i^2I_3\in\mathbb{R}^{3\times 3},
    \end{align}
    
.. math::

    \begin{align}
        R(q)_{ij} = \Big(\sum_{k=j}^N m_k\Big)L_iL_j\hat{q}_i^T\hat{q}_j\in\mathbb{R}^{3\times 3} = R(q)_{ji}^T,\; i<j,
    \end{align}

and 

.. math::

    \begin{align}
        M_{ij} =\Big(\sum_{k={\text{max}}\{i,j\}}^N m_k\Big)L_iL_j I_3\in\mathbb{R}^{3\times 3}.
    \end{align}
    
Equations :ref:`(1) <eq:1>` and :ref:`(2) <eq:2>` define the dynamics of the N-fold pendulum, and hence a vector field :math:`F\in\mathfrak{X}((TS^2)^N)`. We now find a function :math:`f:(TS^2)^N\rightarrow \mathfrak{se}(3)^N` such that

.. math::

    \begin{align}
        \Psi_*(f(m))\vert_m = F\vert_m,\;\;\forall m\in (TS^2)^N.
    \end{align}

Since :math:`R(q)` defines a linear invertible map (see `(Celledoni, Çokaj, Leone, Murari and Owren, (2021) International Journal of Computer Mathematics) <https://doi.org/10.1080/00207160.2021.1966772>`_).

.. math::

    \begin{align}
        A_{q}:T_{q_1}S^2\times ... \times T_{q_N}S^2 \rightarrow T_{q_1}S^2 \times ... \times T_{q_N}S^2,\quad A_q(\omega):=R(q)\omega,
    \end{align}
    
we can rewrite the ODEs for the angular velocities as follows:

.. math::
    :name: eq:3
    
    \begin{align}
        \dot{\omega}= A_{q}^{-1}\left(\begin{bmatrix}r_1\\ \vdots \\ r_N \end{bmatrix}\right) =\begin{bmatrix} h_1(q,\omega) \\ \vdots \\ h_N(q,\omega)\end{bmatrix} = \begin{bmatrix} a_1(q,\omega)\times q_1 \\ \vdots \\ a_N(q,\omega)\times q_N \end{bmatrix}.
    \end{align}
   
In equation :ref:`(3) <eq:3>` the :math:`r_i-s` are defined as in :ref:`(2) <eq:2>` ,
and :math:`a_1,...,a_N:(TS^2)^N\rightarrow \mathbb{R}^3` can be defined as :math:`a_i(q,\omega):=q_i\times h_i(q,\omega)`. Thus, the map :math:`f` is given by

.. math::

    \begin{align}
        f(q,\omega) = \begin{bmatrix}
        \omega_1 \\
        q_1\times h_1(q,\omega) \\ \vdots \\ \omega_N \\ q_N\times h_N(q,\omega)
        \end{bmatrix}\in\mathfrak{se}(3)^N\simeq \mathbb{R}^{6N}.
   \end{align}
