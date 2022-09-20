
.. _lie_group_integrators:

=====================
Integration methods
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
} such that

.. math::
    :name: eq:

    \begin{align}
        F|_m = \left.\Psi_*(f(m))\right|_m,\quad\text{for all}\; m\in \mathcal{M}
    \end{align}

This is the original tool from `(H. Munthe-Kaas, (1999) Appl. Numer. Math 29) <https://doi.org/10.1016/S0168-9274(98)00030-0>`_ for representing a vector field on a manifold with a group action.

.. _two_classes_lie_groups:

Two classes of Lie group integrators
====================================

The simplest numerical integrator for linear spaces is the explicit Euler method. 
Given an initial value problem :math:`\dot{y}=F(y)`, :math:`y(0)=y_0` the method is defined as :math:`y_{n+1}=y_n + hF(y_n)` for some stepsize :math:`h`. 
In the spirit of the previous section, one could think of
the Euler method as the :math:`h` -flow of the constant vector field :math:`F_{y_n}(y)=F(y_n)`, that is

.. math::
    :name: eq:
    
    \begin{align}
        y_{n+1} = \exp(hF_{y_n})\,y_n
    \end{align}

This definition of the Euler method makes sense also when :math:`F` is replaced by a vector field on some manifold. 
In this general situation it is known as the Lie-Euler method.

We briefly consider here two classes of methods known as Runge-Kutta-Munthe-Kaas (RKMK) methods and Commutator-free Lie group methods. 
For a more detailed discussion, the reader can see `(Celledoni, Çokaj, Leone, Murari and Owren, (2021) International Journal of Computer Mathematics) <https://doi.org/10.1080/00207160.2021.1966772>`_ and the references therein.

.. _RKMK:

Runge-Kutta-Munthe-Kaas (RKMK) methods
--------------------------------------

For RKMK methods the underlying idea is to transform the problem from the manifold :math:`\mathcal{M}` to the Lie algebra :math:`\mathfrak{g}`, take a time step, and map the result back to :math:`\mathcal{M}`. The transformation we use is

.. math::
    :name: eq:
    
    \begin{align}
        y(t) = \exp(\sigma(t))\cdot y_0,\quad\sigma(0)=0.
    \end{align}
    
The transformed differential equation for :math:`\sigma(t)` makes use of the derivative of the exponential mapping,

.. math::
    :name: eq:3

    \begin{align}
        \dot{\sigma}(t) = \textrm{dexp}_{\sigma(t)}^{-1} (f(\exp(\sigma(t))\cdot y_0))
    \end{align}

The map :math:`v\mapsto\textrm{dexp}_u(v)` is linear and invertible when :math:`u` belongs to some sufficiently small neighborhood of :math:`0\in\mathfrak{g}`. 
It has an expansion in nested Lie brackets. Using the operator :math:`\textrm{ad}_u(v)=[u,v]` and its powers
:math:`\textrm{ad}_u^2 v=[u,[u,v]]` etc, one can write

.. math::
    :name: eq:4
    
    \begin{align}
        \textrm{dexp}_u(v) = \left.\frac{e^z-1}{z}\right|_{z=\textrm{ad}_u}(v) = v + \frac12[u,v] + \frac16[u,[u,v]] + \cdots
    \end{align}
    
and the inverse is

.. math::
    :name: eq:5
    
    \begin{align}
        \textrm{dexp}_u^{-1}(v) =\left.\frac{z}{e^z-1}\right|_{z=\textrm{ad}_u}(v)= v -\frac12[u,v] + \frac1{12}[u,[u,v]]+\cdots
    \end{align}

The RKMK methods are now obtained simply by applying some standard Runge-Kutta method to the transformed equation :ref:`(3) <eq:3>` with a time step :math:`h`, using initial value :math:`\sigma(0)=0`. This leads to an output :math:`\sigma_1\in\mathfrak{g}` and one simply sets :math:`y_1=\exp(\sigma_1)\cdot y_0`. Then one repeats the procedure replacing :math:`y_0` by :math:`y_1` in the next step etc. While solving :ref:`(3) <eq:3>` one needs to evaluate :math:`\textrm{dexp}_u^{-1}(v)` as a part of the process. This can be done by truncating the series :ref:`(5) <eq:5>` since :math:`\sigma(0)=0` implies that we always evaluate :math:`\textrm{dexp}_u^{-1}` with :math:`u=\mathcal{O}(h)`, and thus, the :math:`k-th` iterated commutator :math:`\textrm{ad}_u^k=\mathcal{O}(h^k)`.
For a given Runge-Kutta method, there are some clever tricks that can be done to minimise the total number of commutators to be included from the expansion of :math:`\textrm{dexp}_u^{-1}v`. A concrete example of an RKMK method is:



.. math::
    
    \begin{align}
       f_{n,1} &= h f(y_n),\\
       f_{n,2} &= h f(\exp(\tfrac{1}{2}f_{n,1}) \cdot y_n), \\
       f_{n,3} &= h f(\exp(\tfrac{1}{2}f_{n,2}-\tfrac{1}{8}[f_{n,1},f_{n,2}])\cdot y_n), \\
       f_{n,4} &= h f(\exp(f_{n,3})\cdot y_n), &  \\
       y_{n+1} &= \exp(\tfrac{1}{6}(f_{n,1}+2f_{n,2}+2f_{n,3}+f_{n,4}-\tfrac12[f_{n,1},f_{n,4}]))\cdot y_n.
    \end{align}
    

The other option is to compute the exact expression for :math:`\textrm{dexp}_u^{-1}(v)` for the particular Lie algebra we use. An expression for  :math:`\textrm{dexp}_u^{-1}(v)` for the Lie algebra :math:`\mathfrak{so}(3)` was shown in `(Celledoni and Owren, (2003) Computer Methods in Applied Mechanics and Engineering) <https://doi.org/10.1016/S0045-7825(02)00520-0>`_. We provide an exact expression for :math:`\textrm{dexp}_u^{-1}(v)` in :math:`\mathfrak{se}(3)` in `(Celledoni, Çokaj, Leone, Murari and Owren, (2021) International Journal of Computer Mathematics) <https://doi.org/10.1080/00207160.2021.1966772>`_

.. _CFmethods:

Commutator-free methods
-----------------------

The second class of Lie group integrators we consider here are the commutator-free methods. As the name suggests, in contrast to RKMK schemes which usually include commutators in the method format, this class of methods does not. These schemes include the Crouch-Grossman methods and they have the format

.. math::

    \begin{align}
        Y_{n,r} &= \exp\left(h\sum_{k}\alpha_{r,J}^k f_{n,k}\right)\cdots \exp\left(h\sum_{k}\alpha_{r,1}^k f_{n,k}\right)\cdot y_n \\
        f_{n,r} &= f(Y_{n,r}) \\[1mm]
        y_{n+1} &= \exp\left(h\sum_k \beta_J^k f_{n,k}\right)\cdots \exp\left(h\sum_k \beta_1^k f_{n,k}\right)\cdot y_n
    \end{align}
    
Here the Runge-Kutta coefficients :math:`\alpha_{r,j}^k`, :math:`\beta_{j}^r` are related to a classical Runge-Kutta scheme with coefficients :math:`a_r^k`, :math:`b_r` in that :math:`a_r^k=\sum_j \alpha_{r,j}^k` and :math:`b_r=\sum_j \beta_{j}^r`. The :math:`\alpha_{r,j}^k`, :math:`\beta_{j}^r` are usually chosen to obtain computationally inexpensive schemes with the highest possible order of convergence. The computational complexity of the above schemes depends on the cost of computing an exponential as well as of evaluating the vector field. Therefore it makes sense to keep the number of exponentials :math:`J` in each stage as low as possible, and possibly also the number of stages :math:`s`. A trick is to select coefficients that make it possible to reuse exponentials from one stage to another. This is perhaps best illustrated through the following example, a generalisation of the classical 4th order Runge-Kutta method.

.. math::
    
    \begin{align}
        \begin{split} \label{cf4}
        Y_{n,1} &= y_n \\
        Y_{n,2} &=\exp(\tfrac12 hf_{n,1}) \cdot y_n \\
        Y_{n,3} &= \exp(\tfrac12 hf_{n,2}) \cdot y_n \\
        Y_{n,4} &= \exp(h f_{n,3}-\tfrac12 h f_{n,1}) \cdot Y_{n,2} \\
        y_{n+\frac12} &=\exp(\tfrac{1}{12}h(3f_{n,1}+2f_{n,2}+2f_{n,3}-f_{n,4})) \cdot y_n \\
        y_{n+1} &=\exp(\tfrac{1}{12}h(-f_{n,1}+2f_{n,2}+2f_{n,3}+ 3f_{n,4})) \cdot y_{n+\frac12}
        \end{split}
    \end{align}

where :math:`f_{n,i}=f(Y_{n,i})`. Here, we see that one exponential is saved in computing :math:`Y_{n,4}` by making use of :math:`Y_{n,2}`.

In our `code <https://github.com/THREAD-3-2/lie_group_integrators_double_spherical_pendulum/tree/main/src>`_ we have tested both types of methods discussed here. 
We have performed our numerical experiments on the example of `the double spherical pendulum <https://thread-3-2.github.io/lie_group_integrators_double_spherical_pendulum/double_sph_pend.html>`_.
Experiments show that Lie group integrators allow to keep the evolution of the solution in the correct manifold.
We show the convergence rate of all the Lie group integrators tested on this model and we check how they behave in terms of preserving the configuration manifold and the phase space.
The analysis is completed with a comparison with the classical Runge-Kutta 4 and with ODE45 of MATLAB. 
The Lie group integrators used to obtain the experiments are Lie-Euler, Lie-Euler-Heun, three versions of Runge-Kutta-Munthe-Kaas methods of order four and one of order three. 
Their implementation can be found `here <https://github.com/THREAD-3-2/lie_group_integrators_double_spherical_pendulum/tree/main/src/integrators>`_.
Unlike classical numerical integrators like the one implemented in ODE45 or the Runge-Kutta 4, the Lie group methods preserve the configuration manifold and the phase space to a high accuracy. 
 
