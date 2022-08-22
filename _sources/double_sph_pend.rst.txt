.. _double_sph_pend:

=============================
The double spherical pendulum
=============================

In the paper `(Celledoni, Çokaj, Leone, Murari and Owren, (2021) International Journal of Computer Mathematics) <https://doi.org/10.1080/00207160.2021.1966772>`_ we discuss in detail the N-fold 3D pendulum. 
We consider here the special case of the double spherical pendulum,  :math:`N = 2`. 
We then test our Lie group integrators on this particular case, showing the rate of convergence, the preservation of the geometry of the manifold :math:`S^2` and the configuration space :math:`T_{q_{i}(t)}S^2`.


.. _trans_action:

Transitive group action on (TS^2)^2
-----------------------------------

We characterize a transitive action for :math:`(TS^2)^N`, starting with the case :math:`N=1` and generalizing it to :math:`N>1` . The action we consider is based on the identification between :math:`\mathfrak{se}(3)`, the Lie algebra of :math:`SE(3)`, and :math:`\mathbb{R}^6`. We start from the Ad-action of :math:`SE(3)` on :math:`\mathfrak{se}(3)` (see \cite{Holm}), which writes

.. math::
    :name: eq: 
    
    \begin{align}
\Ad : SE(3)\times \mathfrak{se}(3) \rightarrow \mathfrak{se}(3),
\end{align}

.. math::
    :name: eq: 
    
    \begin{align}
\Ad((R,r),(u,v)) = (Ru,Rv+\hat{r}Ru).
\end{align}
Since :math:`\mathfrak{se}(3)\simeq \mathbb{R}^6`, the Ad-action allows us to define the following Lie group action on :math:`\mathbb{R}^6`

.. math::
    :name: eq: 
    
    \begin{align}
\psi: SE(3)\times\mathbb{R}^6\rightarrow \mathbb{R}^6,\;\;\psi((R,r),(u,v)) = (Ru,Rv+\hat{r}Ru).
\end{align}
We can think of :math:`\psi` as a Lie group action on :math:`TS^2` since, for any :math:`q\in\mathbb{R}^3`, it maps points of

.. math::
    :name: eq: 
    
    \begin{align}TS_{|q|}^2:=\{(\Tilde{q},\Tilde{\omega})\in \mathbb{R}^3\times\mathbb{R}^3:\; \Tilde{\omega}^T\Tilde{q}=0,\;|\Tilde{q}|=|q|\}\subset \mathbb{R}^6
\end{align}
into other points of :math:`TS_{|q|}^2`. Moreover, with standard arguments (see \cite{olver2000applications}), it is possible to prove that the orbit of a generic point :math:`m=(q,\omega)\in\mathbb{R}^6` with :math:`\omega^Tq=0` coincides with

.. math::
    :name: eq: 
    
    \begin{align}
\RE{\text{Orb}}(m)=TS_{|q|}^2.
\end{align}
In particular, when :math:`q\in\mathbb{R}^3` is a unit vector (i.e. :math:`q\in S^2`), :math:`\psi` allows us to define a transitive Lie group action on :math:`TS^2=TS_{|q|=1}^2` which writes

.. math::
    :name: eq: 
    
    \begin{align}
\psi : SE(3)\times TS^2 \rightarrow TS^2
\end{align}

.. math::
    :name: eq: 
    
    \begin{align}
\psi((A,a),(q,\omega)) := \psi_{(A,a)}(q,\omega) =  (Aq,A\omega + \hat{a}Aq)=(\bar{q},\bar{\omega}).
\end{align}
To conclude the description of the action, we report here its infinitesimal generator which is fundamental in the Lie group integrators setting

.. math::
    :name: eq: 
    
    \begin{align}
\left.\infgen((u,v))\right|_{(q,\omega)} =(\hat{u}q,\hat{u}\omega + \hat{v}q).
\end{align}
We can extend this construction to the case :math:`N>1` in a natural way, i.e. through the action of a Lie group obtained from cartesian products of :math:`SE(3)` and equipped with the direct product structure. More precisely, we consider the group :math:`G=(SE(3))^N` and by direct product structure we mean that for any pair of elements 

.. math::
    :name: eq: 
    
    \begin{align}
    \delta^{(1)}=(\delta^{(1)}_1,...,\delta^{(1)}_N),\quad \delta^{(2)}=(\delta^{(2)}_1,...,\delta^{(2)}_N)\in G,
    \end{align}
    
    denoted with :math:`*` the semidirect product of :math:`SE(3)`, we define the product :math:`\circ` on :math:`G` as

.. math::
    :name: eq: 
    
    \begin{align}
\delta^{(1)}\circ \delta^{(2)} := (\delta^{(1)}_1 * \delta^{(2)}_1,...,\delta^{(1)}_N * \delta^{(2)}_N)\in G.
\end{align}
With this group structure defined, we can generalize the action introduced for :math:`N=1` to larger :math:`N`s as follows

.. math::
    :name: eq: 
    
    \begin{align}
\psi : (SE(3))^N\times (TS^2)^N \rightarrow (TS^2)^N,
\end{align}

.. math::
    :name: eq: 
    
    \begin{align}
\begin{split}
\psi&((A_1,a_1,...,A_N,a_n),(q_1,\omega_1,...,q_N,\omega_N)) =\\ &=(A_1q_1,A_1\omega_1+\hat{a}_1A_1q_1,...,A_Nq_N,A_N\omega_N+\hat{a}_NA_Nq_N),
\end{split}
\end{align}
whose infinitesimal generator writes

.. math::
    :name: eq: 
    
    \begin{align}
\infgen(\xi)\vert_m =(\hat{u}_1q_1,\hat{u}_1\omega_1+\hat{v}_1q_1,...,\hat{u}_Nq_N,\hat{u}_N\omega_N+\hat{v}_Nq_N),
\end{align}
where :math:`\xi=[u_1,v_1,...,u_N,v_N]\in\mathfrak{se}(3)^N` and :math:`m=(q_1,\omega_1,...,q_N,\omega_N)\in (TS^2)^N`.
We have now the only group action we need to deal with the :math:`N-`fold spherical pendulum. In the following part of this section we work on the vector field describing the dynamics and adapt it to the Lie group integrators setting.

\subsection{Full chain}
We consider the vector field :math:`F\in\mathfrak{X}((TS^2)^N)`, describing the dynamics of the :math:`N`-fold 3D pendulum, and we express it in terms of the infinitesimal generator of the action defined above. More precisely, we find a function :math:`F:(TS^2)^N\rightarrow \mathfrak{se}(3)^N` such that

.. math::
    :name: eq: 
    
    \begin{align}
\infgen(f(m))\vert_m = F\vert_m,\;\;\forall m\in (TS^2)^N.
\end{align}
We omit the derivation of :math:`F` starting from the Lagrangian of the system, which can be found in the section devoted to mechanical systems on :math:`(S^2)^N` of \cite{lee18gfo}. 
%\colorbox{BurntOrange}{Davide: What if we write down the Lagrangian only, then}
The configuration manifold of the system is :math:`(S^2)^N`, while the Lagrangian, expressed in terms of the variables :math:`(q_1,\omega_1,...,q_N,\omega_N)\in (TS^2)^N`, writes

.. math::
    :name: eq: 
    
    \begin{align}
L(q,\omega) = T(q,\omega)-U(q) =\frac{1}{2}\sum_{i,j=1}^N\Big(M_{ij}\omega_i^T\hat{q}_i^T\hat{q}_j\omega_j\Big) - \sum_{i=1}^N\Big(\sum_{j=i}^N m_j\Big)gL_ie_3^Tq_i,
\end{align}
where

.. math::
    :name: eq: 
    
    \begin{align}
M_{ij} =\Big(\sum_{k=\RE{\text{max}}\{i,j\}}^N m_k\Big)L_iL_j I_3\in\mathbb{R}^{3\times 3}
\end{align}
is the inertia matrix of the system\RE{, :math:`I_3` is the :math:`3\times 3` identity matrix,} and :math:`e_3 = [0,0,1]^T`. Noticing that when :math:`i=j` we get

.. math::
    :name: eq: 
    
    \begin{align}
\omega_i^T\hat{q}_i^T\hat{q}_i\omega_i = \omega_i^T(I_3-q_iq_i^T)\omega_i = \omega_i^T\omega_i,
\end{align}
we simplify the notation writing 

.. math::
    :name: eq: 
    
    \begin{align}
T(q,\omega) = \frac{1}{2}\sum_{i,j=1}^N\Big(\omega_i^TR(q)_{ij}\omega_j\Big)
\end{align}
where :math:`R(q)\in\mathbb{R}^{3N\times 3N}` is a symmetric block matrix defined as

.. math::
    :name: eq: 
    
    \begin{align}
R(q)_{ii} = \Big(\sum_{j=i}^Nm_j\Big)L_i^2I_3\in\mathbb{R}^{3\times 3},
\end{align}

.. math::
    :name: eq: 
    
    \begin{align}
R(q)_{ij} = \Big(\sum_{k=j}^N m_k\Big)L_iL_j\hat{q}_i^T\hat{q}_j\in\mathbb{R}^{3\times 3} = R(q)_{ji}^T,\; i<j.
\end{align}


.. _eom_vec_field:

The equations of motion and the vector field
--------------------------------------------

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

This map guarantees that if we rewrite the pair of equations for the angular velocities in :ref:`(2) <eq:rq>` as

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