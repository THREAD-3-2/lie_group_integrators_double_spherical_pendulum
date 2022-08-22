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