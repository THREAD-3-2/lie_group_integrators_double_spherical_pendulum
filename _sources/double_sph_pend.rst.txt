.. _double_sph_pend:

=============================
The double spherical pendulum
=============================

In the paper `(Celledoni, Çokaj, Leone, Murari and Owren, (2021) International Journal of Computer Mathematics) <https://doi.org/10.1080/00207160.2021.1966772>`_ we discuss in detail the N-fold spherical pendulum, which is an example of the initial value problem discussed in the first `section <https://thread-3-2.github.io/lie_group_integrators_double_spherical_pendulum/lie_group_integrators.html>`_. 
We consider here the special case of the double spherical pendulum,  :math:`N = 2`. 
We then test our Lie group integrators on this particular case, showing the rate of convergence, the preservation of the geometry of the manifold :math:`S^2` and the phase space :math:`T_{q_{i}(t)}S^2`.

The Lagrangian we consider is a function from :math:`(TS^2)^2` to :math:`\mathbb{R}`. Instead of the coordinates :math:`(q_1, q_2,\dot{q}_1, \dot{q}_2)`, where :math:`\dot{q}_i\in T_{q_i}S^2`, we choose to work with the angular velocities. 
Precisely, 

.. math::
    :name: eq: 
    
    \begin{align}
        T_{q_i}S^2 = \{v\in\mathbb{R}^3:\;v^Tq_i=0\} = \langle q_i\rangle ^{\perp} \subset \mathbb{R}^3, \quad i = 1,2
    \end{align}

and hence for any :math:`\dot{q}_i\in T_{q_i}S^2` there exist :math:`\omega_i\in\mathbb{R}^3` such that :math:`\dot{q}_i=\omega_i\times q_i`, which can be interpreted as the angular velocity of :math:`q_i`. 
So we can assume without loss of generality that :math:`\omega_i^Tq_i=0` (i.e. :math:`\omega_i\in T_{q_i}S^2`) and pass to the coordinates :math:`(q_1,\omega_1,q_2,\omega_2)\in (TS^2)^2` to describe the dynamics.  
We denote with :math:`m_1, m_2` the masses of the pendulums and with :math:`L_1, L_2` their lengths.


.. _trans_action:

Transitive group action on :math:`(TS^2)^2`
-------------------------------------------

In `(Celledoni, Çokaj, Leone, Murari and Owren, (2021) International Journal of Computer Mathematics) <https://doi.org/10.1080/00207160.2021.1966772>`_ we characterize a transitive action for :math:`(TS^2)^N`, starting with the case :math:`N=1` and generalizing it to :math:`N>1` . 
The action we consider is based on the identification between :math:`\mathfrak{se}(3)`, the Lie algebra of :math:`SE(3)`, and :math:`\mathbb{R}^6`. We start from the Ad-action of :math:`SE(3)` on :math:`\mathfrak{se}(3)`, which is

.. math::
    :name: eq: 
    
    \begin{align}
        \textrm{Ad} : SE(3)\times \mathfrak{se}(3) \rightarrow \mathfrak{se}(3),
    \end{align}

.. math::
    :name: eq: 
    
    \begin{align}
        \textrm{Ad}((R,r),(u,v)) = (Ru,Rv+\hat{r}Ru).
    \end{align}

Since :math:`\mathfrak{se}(3)\simeq \mathbb{R}^6`, the Ad-action allows us to define the following Lie group action on :math:`\mathbb{R}^6`

.. math::
    :name: eq: 
    
    \begin{align}
        \Psi: SE(3)\times\mathbb{R}^6\rightarrow \mathbb{R}^6,\;\;\Psi((R,r),(u,v)) = (Ru,Rv+\hat{r}Ru).
    \end{align}

We can think of :math:`\Psi` as a Lie group action on :math:`TS^2` since, for any :math:`q\in\mathbb{R}^3`, it maps points of

.. math::
    :name: eq: 
    
    \begin{align}
        TS_{|q|}^2:=\{(\tilde{q},\tilde{\omega})\in \mathbb{R}^3\times\mathbb{R}^3:\; \tilde{\omega}^T\tilde{q}=0,\;|\tilde{q}|=|q|\}\subset \mathbb{R}^6
    \end{align}

into other points of :math:`TS_{|q|}^2`. In particular, when :math:`q\in\mathbb{R}^3` is a unit vector (i.e. :math:`q\in S^2`), :math:`\Psi` allows us to define a transitive Lie group action on :math:`TS^2=TS_{|q|=1}^2` which is

.. math::
    :name: eq: 
    
    \begin{align}
        \Psi : SE(3)\times TS^2 \rightarrow TS^2
    \end{align}

.. math::
    :name: eq: 
    
    \begin{align}
        \Psi((A,a),(q,\omega)) := \Psi_{(A,a)}(q,\omega) =  (Aq,A\omega + \hat{a}Aq)=(\bar{q},\bar{\omega}).
    \end{align}

To conclude the description of the action, we report here its infinitesimal generator which is fundamental in the Lie group integrators setting

.. math::
    :name: eq: 
    
    \begin{align}
        \Psi_*((u,v))\|_{(q,\omega)} =(\hat{u}q,\hat{u}\omega + \hat{v}q).
    \end{align}

We can extend this construction to the case :math:`N>1` in a natural way, i.e. through the action of a Lie group obtained from cartesian products of :math:`SE(3)` and equipped with the direct product structure. 

Here we limit ourselves to the case :math:`N=2` for which we also show numerical experiments. 
More precisely, we consider the group :math:`G=(SE(3))^2` and by direct product structure we mean that for any pair of elements 

.. math::
    :name: eq: 
    
    \begin{align}
        \delta^{(1)}=(\delta^{(1)}_1, \delta^{(1)}_2), \delta^{(2)}=(\delta^{(2)}_1, \delta^{(2)}_2)\in G,
    \end{align}
    
denoted with :math:`*` the semidirect product of :math:`SE(3)`, we define the product :math:`\circ` on :math:`G` as

.. math::
    :name: eq: 
    
    \begin{align}
        \delta^{(1)}\circ \delta^{(2)} := (\delta^{(1)}_1 * \delta^{(2)}_1, \delta^{(1)}_2 * \delta^{(2)}_2)\in G.
    \end{align}

With this group structure defined, we can write the action follows

.. math::
    :name: eq: 
    
    \begin{align}
        \Psi : (SE(3))^2\times (TS^2)^2 \rightarrow (TS^2)^2,
    \end{align}

.. math::
    :name: eq: 
    
    \begin{align}
        \begin{split}
        \Psi&((A_1,a_1, A_2,a_2),(q_1,\omega_1, q_2,\omega_2)) =\\ &=(A_1q_1,A_1\omega_1+\hat{a}_1A_1q_1, A_2q_2,A_2\omega_2+\hat{a}_2A_2q_2),
        \end{split}
    \end{align}

whose infinitesimal generator is

.. math::
    :name: eq: 
    
    \begin{align}
        \Psi_*(\xi)\vert_m =(\hat{u}_1q_1,\hat{u}_1\omega_1+\hat{v}_1q_1, \hat{u}_2q_2,\hat{u}_2\omega_2+\hat{v}_2q_2),
    \end{align}

where :math:`\xi=[u_1,v_1, u_2,v_2]\in\mathfrak{se}(3)^2` and :math:`m=(q_1,\omega_1, q_2,\omega_2)\in (TS^2)^2`.
We have now the only group action we need to deal with the double spherical pendulum. In the following part of this section we work on the vector field describing the dynamics and adapt it to the Lie group integrators setting.


The equations of motion and the vector field
--------------------------------------------

We consider the vector field :math:`F\in\mathfrak{X}((TS^2)^2)`, describing the dynamics of the double spherical pendulum, and we express it in terms of the infinitesimal generator of the action defined above. 
More precisely, we find a function :math:`F:(TS^2)^2\rightarrow \mathfrak{se}(3)^2` such that

.. math::
    :name: eq: 
    
    \begin{align}
        \Psi_*(f(m))\vert_m = F\vert_m,\;\;\forall m\in (TS^2)^2.
    \end{align}

The derivation of :math:`F` starting from the Lagrangian of the system can be found in the section devoted to mechanical systems on :math:`(S^2)^2` of `(Lee, Leok and McClamroch, (2018)) <https://doi.org/10.1007/978-3-319-56953-6>`_. 
The configuration manifold of the system is :math:`(S^2)^2`, while the Lagrangian, expressed in terms of the variables :math:`(q_1,\omega_1, q_2,\omega_2)\in (TS^2)^2`, is

.. math::
    :name: eq: 
    
    \begin{align}
        L(q,\omega) = T(q,\omega)-U(q) =\frac{1}{2}\sum_{i,j=1}^2\Big(M_{ij}\omega_i^T\hat{q}_i^T\hat{q}_j\omega_j\Big) - \sum_{i=1}^2\Big(\sum_{j=i}^2 m_j\Big)gL_ie_3^Tq_i,
    \end{align}

where

.. math::
    :name: eq: 
    
    \begin{align}
        M_{ij} =\Big(\sum_{k=\textrm{max}\{i,j\}}^2 m_k\Big)L_iL_j I_3\in\mathbb{R}^{3\times 3}
    \end{align}

is the inertia matrix of the system, :math:`I_3` is the :math:`3\times 3` identity matrix, and :math:`e_3 = [0,0,1]^T`. Noticing that when :math:`i=j` we get

.. math::
    :name: eq: 
    
    \begin{align}
        \omega_i^T\hat{q}_i^T\hat{q}_i\omega_i = \omega_i^T(I_3-q_iq_i^T)\omega_i = \omega_i^T\omega_i,
    \end{align}

we simplify the notation writing 

.. math::
    :name: eq: 
    
    \begin{align}
        T(q,\omega) = \frac{1}{2}\sum_{i,j=1}^2\Big(\omega_i^TR(q)_{ij}\omega_j\Big)
    \end{align}

where :math:`R(q)\in\mathbb{R}^{6\times 6}` is a symmetric block matrix defined as

.. math::
    :name: eq: 
    
    \begin{align}
        R(q)_{ii} = \Big(\sum_{j=i}^2m_j\Big)L_i^2I_3\in\mathbb{R}^{3\times 3},
    \end{align}


.. math::
    :name: eq: 
    
    \begin{align}
        R(q)_{ij} = \Big(\sum_{k=j}^2 m_k\Big)L_iL_j\hat{q}_i^T\hat{q}_j\in\mathbb{R}^{3\times 3} = R(q)_{ji}^T,\; i<j.
    \end{align}


Precisely, the equations of motion write:

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
    :name: eq:22
    
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

This map guarantees that if we rewrite the pair of equations for the angular velocities in :ref:`(22) <eq:22>` as

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
        \Psi_*(f(q,\omega))\vert_{(q,\omega)}=F\vert_{(q,\omega)}
    \end{align}

through the function

.. math::
    :name: eq: 
    
    \begin{align}
        f : TS^2\times TS^2\rightarrow \mathfrak{se}(3)\times\mathfrak{se}(3)\simeq \mathbb{R}^{12},\;\;(q,\omega)\rightarrow (\omega_1, q_1\times h_1, \omega_2,q_2\times h_2).
    \end{align}
