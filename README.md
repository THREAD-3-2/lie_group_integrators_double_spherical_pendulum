# Lie group integrators double spherical pendulum

![CI](https://github.com/THREAD-3-2/lie_group_integrators_double_spherical_pendulum/workflows/CI/badge.svg)
[![documentation](https://img.shields.io/badge/docs-passing-<COLOR>.svg)](https://THREAD-3-2.github.io/lie_group_integrators_double_spherical_pendulum/)

Different ways of applying Lie group integrators to simulating the dynamics of mechanical multi-body systems are discussed here.
The point of departure is the formulation of the models as differential equations on manifolds. Assuming to be given either a Lie group acting transitively on the manifold $\mathcal{M}$ or a set of frame vector fields on $\mathcal{M}$, we use them to describe the mechanical system and further to build the numerical integrator. We consider schemes of the types commonly known as Runge-Kutta-Munthe-Kaas methods and Commutator-free Lie group methods. 
Lie group integrators applied to the example of the double spherical pendulum. The numerical experiments are performed on the example of the double spherical pendulum. We show the convergence rate of all the Lie group integrators tested on this model and we check how they behave in terms of preserving the configuration manifold and the phase space. The analysis is completed with a comparison with the classical Runge–Kutta 4 and with ODE45 of MATLAB. The Lie group integrators used to obtain the experiments are Lie Euler, Lie Euler Heun, three versions of Runge–Kutta–Munthe–Kaas methods of order four and one of order three. Unlike classical numerical integrators like the one implemented in ODE45 or the Runge–Kutta 4, the Lie group methods preserve the configuration manifold and the phase space to a high accuracy.

Matlab code.
