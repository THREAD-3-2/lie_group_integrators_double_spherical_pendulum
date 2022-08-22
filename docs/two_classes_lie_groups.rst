.. _two_classes_lie_groups:

====================================
Two classes of Lie group integrators
====================================

We consider here two classes of methods known as Runge--Kutta--Munthe--Kaas (RKMK) methods and Commutator-free Lie group methods. 

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
    :name: eq:
    
    \begin{align}
        \textrm{dexp}_u(v) = \left.\frac{e^z-1}{z}\right|_{z=\textrm{ad}_u}(v) = v + \frac12[u,v] + \frac16[u,[u,v]] + \cdots
    \end{align}
    
and the inverse is

.. math::
    :name: eq:2
    
    \begin{align}
        \textrm{dexp}_u^{-1}(v) =\left.\frac{z}{e^z-1}\right|_{z=\textrm{ad}_u}(v)= v -\frac12[u,v] + \frac1{12}[u,[u,v]]+\cdots
    \end{align}

The RKMK methods are now obtained simply by applying some standard Runge--Kutta method to the transformed equation :ref:`(1) <eq:1>` with a time step :math:`h`, using initial value :math:`\sigma(0)=0`. This leads to an output :math:`\sigma_1\in\mathfrak{g}` and one simply sets :math:`y_1=\exp(\sigma_1)\cdot y_0`. Then one repeats the procedure replacing :math:`y_0` by :math:`y_1` in the next step etc. While solving :ref:`(1) <eq:1>` one needs to evaluate :math:`\textrm{dexp}_u^{-1}(v)` as a part of the process. This can be done by truncating the series :ref:`(2) <eq:2>` since :math:`\sigma(0)=0` implies that we always evaluate :math:`\textrm{dexp}_u^{-1}` with :math:`u=\mathcal{O}(h)`, and thus, the :math:`k-th` iterated commutator :math:`\textrm{ad}_u^k=\mathcal{O}(h^k)`.
For a given Runge--Kutta method, there are some clever tricks that can be done to minimise the total number of commutators to be included from the expansion of :math:`\textrm{dexp}_u^{-1}v`. We give here one concrete example of an RKMK method:



.. math::
    
    \begin{align}
       f_{n,1} &= h f(y_n),\\
       f_{n,2} &= h f(\exp(\tfrac{1}{2}f_{n,1}) \cdot y_n), \\
       f_{n,3} &= h f(\exp(\tfrac{1}{2}f_{n,2}-\tfrac{1}{8}[f_{n,1},f_{n,2}])\cdot y_n), \\
       f_{n,4} &= h f(\exp(f_{n,3})\cdot y_n), &  \\
       y_{n+1} &= \exp(\tfrac{1}{6}(f_{n,1}+2f_{n,2}+2f_{n,3}+f_{n,4}-\tfrac12[f_{n,1},f_{n,4}]))\cdot y_n.
    \end{align}



The other option is to compute the exact expression for :math:`\textrm{dexp}_u^{-1}(v)` for the particular Lie algebra we use. For instance, it was shown in \cite{celledoni03lgm} that for the Lie algebra :math:`\mathfrak{so}(3)` one has

.. math::
    
    \begin{align}
        \textrm{dexp}_u^{-1}(v)=v - \frac12 u\times v + \alpha^{-2}(1-\tfrac{\alpha}{2}\cot\tfrac{\alpha}{2})\; u\times (u\times v)
    \end{align}
    
We will present the corresponding formula for :math:`\mathfrak{se}(3)` in Section~\ref{dexpinvse3}.


.. _CFmethods:

Commutator-free methods
-----------------------

The second class of Lie group integrators to be considered here are the commutator-free methods, named this way in \cite{celledoni03cfl} to emphasize the contrast to RKMK schemes which usually include commutators in the method format. These schemes include the Crouch-Grossman methods \cite{crouch93nio} and they have the format

.. math::

    \begin{align}
        Y_{n,r} &= \exp\left(h\sum_{k}\alpha_{r,J}^k f_{n,k}\right)\cdots \exp\left(h\sum_{k}\alpha_{r,1}^k f_{n,k}\right)\cdot y_n \\
        f_{n,r} &= f(Y_{n,r}) \\[1mm]
        y_{n+1} &= \exp\left(h\sum_k \beta_J^k f_{n,k}\right)\cdots \exp\left(h\sum_k \beta_1^k f_{n,k}\right)\cdot y_n
    \end{align}
    
Here the Runge--Kutta coefficients :math:`\alpha_{r,j}^k`, :math:`\beta_{j}^r` are related to a classical Runge--Kutta scheme with coefficients :math:`a_r^k`, :math:`b_r` in that :math:`a_r^k=\sum_j \alpha_{r,j}^k` and :math:`b_r=\sum_j \beta_{j}^r`. The :math:`\alpha_{r,j}^k`, :math:`\beta_{j}^r` are usually chosen to obtain computationally inexpensive schemes with the highest possible order of convergence. The computational complexity of the above schemes depends on the cost of computing an exponential as well as of evaluating the vector field. Therefore it makes sense to keep the number of exponentials :math:`J` in each stage as low as possible, and possibly also the number of stages :math:`s`. A trick proposed in \cite{celledoni03cfl} was to select coefficients \RE{that make it possible to reuse exponentials} from one stage to another. This is perhaps best illustrated through the following example from \cite{celledoni03cfl}, a generalisation of the classical 4th order Runge--Kutta method.

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