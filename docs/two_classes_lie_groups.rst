.. _two_classes_lie_groups:

====================================
Two classes of Lie group integrators
====================================

The simplest numerical integrator for linear spaces is the explicit Euler method. 
Given an initial value problem :math:`\dot{y}=F(y)`, :math:`y(0)=y_0` the method is defined as :math:`y_{n+1}=y_n + hF(y_n)` for some stepsize :math:`h`. 
In the spirit of the previous section, one could think of
the Euler method as the :math:`h` -flow of the constant vector field :math:`F_{y_n}(y)=F(y_n)`, that is
.. math::
    
    \begin{align}
        y_{n+1} = \exp(hF_{y_n})\,y_n
    \end{align}

This definition of the Euler method makes sense also when :math:`F` is replaced by a vector field on some manifold. 
In this general situation it is known as the Lie--Euler method.

We briefly consider here two classes of methods known as Runge--Kutta--Munthe--Kaas (RKMK) methods and Commutator-free Lie group methods. 
For a more detailed discussion, the reader can see `(Celledoni, Çokaj, Leone, Murari and Owren, (2021) International Journal of Computer Mathematics) <https://doi.org/10.1080/00207160.2021.1966772>`_ and the references therein.

.. _RKMK:

Runge-Kutta-Munthe-Kaas (RKMK) methods
--------------------------------------

For RKMK methods the underlying idea is to transform the problem from the manifold :math:`\mathcal{M}` to the Lie algebra :math:`\mathfrak{g}`, take a time step, and map the result back to :math:`\mathcal{M}`. The transformation we use is

.. math::
    
    \begin{align}
        y(t) = \exp(\sigma(t))\cdot y_0,\quad\sigma(0)=0.
    \end{align}
    
The transformed differential equation for :math:`\sigma(t)` makes use of the derivative of the exponential mapping,

.. math::
    :name: eq:1

    \begin{align}
        \dot{\sigma}(t) = \textrm{dexp}_{\sigma(t)}^{-1} (f(\exp(\sigma(t))\cdot y_0))
    \end{align}

The map :math:`v\mapsto\textrm{dexp}_u(v)` is linear and invertible when :math:`u` belongs to some sufficiently small neighborhood of :math:`0\in\mathfrak{g}`. 
It has an expansion in nested Lie brackets. Using the operator :math:`\textrm{ad}_u(v)=[u,v]` and its powers
:math:`\textrm{ad}_u^2 v=[u,[u,v]]` etc, one can write

.. math::
    :name: eq:2
    
    \begin{align}
        \textrm{dexp}_u(v) = \left.\frac{e^z-1}{z}\right|_{z=\textrm{ad}_u}(v) = v + \frac12[u,v] + \frac16[u,[u,v]] + \cdots
    \end{align}
    
and the inverse is

.. math::
    :name: eq:3
    
    \begin{align}
        \textrm{dexp}_u^{-1}(v) =\left.\frac{z}{e^z-1}\right|_{z=\textrm{ad}_u}(v)= v -\frac12[u,v] + \frac1{12}[u,[u,v]]+\cdots
    \end{align}

The RKMK methods are now obtained simply by applying some standard Runge--Kutta method to the transformed equation :ref:`(1) <eq:1>` with a time step :math:`h`, using initial value :math:`\sigma(0)=0`. This leads to an output :math:`\sigma_1\in\mathfrak{g}` and one simply sets :math:`y_1=\exp(\sigma_1)\cdot y_0`. Then one repeats the procedure replacing :math:`y_0` by :math:`y_1` in the next step etc. While solving :ref:`(1) <eq:1>` one needs to evaluate :math:`\textrm{dexp}_u^{-1}(v)` as a part of the process. This can be done by truncating the series :ref:`(3) <eq:3>` since :math:`\sigma(0)=0` implies that we always evaluate :math:`\textrm{dexp}_u^{-1}` with :math:`u=\mathcal{O}(h)`, and thus, the :math:`k-th` iterated commutator :math:`\textrm{ad}_u^k=\mathcal{O}(h^k)`.
For a given Runge--Kutta method, there are some clever tricks that can be done to minimise the total number of commutators to be included from the expansion of :math:`\textrm{dexp}_u^{-1}v`. A concrete example of an RKMK method is:



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
    
Here the Runge--Kutta coefficients :math:`\alpha_{r,j}^k`, :math:`\beta_{j}^r` are related to a classical Runge--Kutta scheme with coefficients :math:`a_r^k`, :math:`b_r` in that :math:`a_r^k=\sum_j \alpha_{r,j}^k` and :math:`b_r=\sum_j \beta_{j}^r`. The :math:`\alpha_{r,j}^k`, :math:`\beta_{j}^r` are usually chosen to obtain computationally inexpensive schemes with the highest possible order of convergence. The computational complexity of the above schemes depends on the cost of computing an exponential as well as of evaluating the vector field. Therefore it makes sense to keep the number of exponentials :math:`J` in each stage as low as possible, and possibly also the number of stages :math:`s`. A trick is to select coefficients that make it possible to reuse exponentials from one stage to another. This is perhaps best illustrated through the following example, a generalisation of the classical 4th order Runge--Kutta method.

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
