.. _double_sph_pend:

=============================
The double spherical pendulum
=============================

In the paper `(Celledoni, Çokaj, Leone, Murari and Owren, (2021) International Journal of Computer Mathematics) <https://doi.org/10.1080/00207160.2021.1966772>`_ we discuss in detail the N-fold 3D pendulum. 
We consider here the special case of the double spherical pendulum,  :math:`N = 2`. 
We then test our Lie group integrators on this particular case, showing the rate of convergence, the preservation of the geometry of the manifold :math:`S^2` and the configuration space :math:`T_{q_{i}(t)}S^2`.


The equations of motion write:

.. math::
    :name: eq: 
    
    \begin{align}
        \dot{q}_1 = \hat{\omega}_1q_1,\quad \dot{q}_2 = \hat{\omega}_2q_2,
    \end{align}

.. math::
    :name: eq: 
    
    \begin{align}
        R(q)\begin{bmatrix}
        \dot{\omega}_1 \\ \dot{\omega}_2
        \end{bmatrix}= 
        \begin{bmatrix}
        (-m_2L_1L_2|\omega_2|^2\hat{q}_2 + (m_1+m_2)gL_1\hat{e}_3)q_1 \\
        (-m_2L_1L_2|\omega_1|^2\hat{q}_1 + m_2gL_2\hat{e}_3)q_2
        \end{bmatrix},
    \end{align}

where 

.. math::
    :name: eq: 
    
    \begin{align}
        R(q) = \begin{bmatrix}
        (m_1+m_2)L_1^2I_3 & m_2L_1L_2\hat{q}_1^T\hat{q}_2 \\
        m_2L_1L_2\hat{q}_2^T\hat{q}_1 & m_2L_2^2I_3
        \end{bmatrix}.
    \end{align}

As presented above, the matrix :math:`R(q)` defines a linear invertible map of the space :math:`T_{q_1}S^2\times T_{q_2}S^2` onto itself:

.. math::
    :name: eq: 
    
    \begin{align}
        A_{(q_1,q_2)}:T_{q_1}S^2\times T_{q_2}S^2\rightarrow T_{q_1}S^2\times T_{q_2}S^2,\;[\omega_1,\omega_2]^T\rightarrow R(q)[\omega_1,\omega_2]^T.
    \end{align}

We can easily see that it is well defined since

.. math::
    :name: eq: rq
    
    \begin{align}
        R(q)\begin{bmatrix}
        \omega_1 \\ \omega_2
        \end{bmatrix} = \begin{bmatrix}
        (m_1+m_2)L_1^2I_3 & m_2L_1L_2\hat{q}_1^T\hat{q}_2 \\
        m_2L_1L_2\hat{q}_2^T\hat{q}_1 & m_2L_2^2I_3
        \end{bmatrix}\begin{bmatrix}
        \hat{v}_1q_1 \\ \hat{v}_2q_2
        \end{bmatrix} = \begin{bmatrix}
        \hat{r}_1q_1\\ \hat{r}_2q_2 
        \end{bmatrix}\in (TS^2)^2
    \end{align}

with 

.. math::
    :name: eq: 
    
    \begin{align}
        r_1(q,\omega):=(m_1+m_2)L_1^2v_1+m_2L_1L_2\hat{q}_2\hat{v}_2q_2,
    \end{align} 

.. math::
    :name: eq: 
    
    \begin{align} 
        r_2(q,\omega):=m_2L_1L_2\hat{q}_1\hat{v}_1q_1+m_2L_2^2v_2. 
    \end{align}

This map guarantees that if we rewrite the pair of equations for the angular velocities in :ref:`() <eq:rq>` as

.. math::
    :name: eq: 
    
    \begin{align}
        \begin{split}
        \dot{\omega}&= R^{-1}(q)\begin{bmatrix}
        (-m_2L_1L_2|\omega_2|^2\hat{q}_2 + (m_1+m_2)gL_1\hat{e}_3)q_1 \\
        (-m_2L_1L_2|\omega_1|^2\hat{q}_1 + m_2gL_2\hat{e}_3)q_2
        \end{bmatrix}=R^{-1}(q)b=\\
        &=A_{(q_1,q_2)}^{-1}(b)=\begin{bmatrix}
        h_1 \\ h_2
        \end{bmatrix}\in T_{q_1}S^2\times T_{q_2}S^2,
        \end{split}
    \end{align}

then we are assured that there exists a pair of functions :math:`a_1,a_2:TS^2\times TS^2\rightarrow\mathbb{R}^3` such that

.. math::
    :name: eq: 
    
    \begin{align}
        \dot{\omega} = \begin{bmatrix}
        a_1(q,\omega)\times q_1 \\ a_2(q,\omega)\times q_2
        \end{bmatrix} = \begin{bmatrix}
        h_1(q) \\ h_2(q)
        \end{bmatrix}.
    \end{align}

Since we want :math:`a_i\times q_i = h_i`, we just impose :math:`a_i=q_i\times h_i` and hence the whole vector field can be rewritten as

.. math::
    :name: eq: 
    
    \begin{align}
        \begin{bmatrix}
        \dot{q}_1 \\ \dot{\omega}_1 \\ \dot{q}_2 \\ \dot{\omega}_2
        \end{bmatrix} = \begin{bmatrix}
        \omega_1 \times q_1 \\ (q_1\times h_1)\times q_1 \\ \omega_2\times q_2 \\ (q_2\times h_2)\times q_2
        \end{bmatrix} = F\vert_{(q,\omega)},
    \end{align}

with :math:`h_i=h_i(q,\omega)` and

.. math::
    :name: eq: 
    
    \begin{align}
        \begin{bmatrix}
        h_1(q,\omega) \\ h_2(q,\omega)
        \end{bmatrix} = R^{-1}(q)\begin{bmatrix}
        (-m_2L_1L_2|\omega_2|^2\hat{q}_2 + (m_1+m_2)gL_1\hat{e}_3)q_1 \\
        (-m_2L_1L_2|\omega_1|^2\hat{q}_1 + m_2gL_2\hat{e}_3)q_2
        \end{bmatrix}.
    \end{align}

Therefore, we can express the whole vector field in terms of the infinitesimal generator of the action of :math:`SE(3)\times SE(3)` as

.. math::
    :name: eq: 
    
    \begin{align}
        \infgen(f(q,\omega))\vert_{(q,\omega)}=F\vert_{(q,\omega)}
    \end{align}

through the function

.. math::
    :name: eq: 
    
    \begin{align}
        f : TS^2\times TS^2\rightarrow \mathfrak{se}(3)\times\mathfrak{se}(3)\simeq \mathbb{R}^{12},\;\;(q,\omega)\rightarrow (\omega_1, q_1\times h_1, \omega_2,q_2\times h_2).
    \end{align}






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
