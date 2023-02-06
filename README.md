# Stochastic FeedBack
<a href = "https://github.com/zaman13/Brownian-dynamics-in-a-time-varying-force-field/tree/main/Codes"> <img src="https://img.shields.io/badge/Language-C++ & Python-orange" alt="alt text"> </a>

The Discrete Nonlinear Schrödinger Equation (DNLS) is a mathematical model that describes the behavior of nonlinear waves in discrete systems, such as lattices, optical fibers, and Bose-Einstein condensates. The equation incorporates both linear and nonlinear terms, and its behavior is influenced by the interaction between waves and their nonlinear interactions. Unlike the linear Schrödinger equation, the DNLS equation predicts the emergence of localized structures such as solitons and breathers, as well as the possibility of chaos. It is widely used in various fields of physics, engineering, and mathematics to study the dynamics of discrete nonlinear systems.

$$
dc_n = i(c_{n+1} + c_{n-1} -2c_n)dt + i\alpha \left | c_n \right |^2c_ndt
$$

# Results


In the simulation, I have utilized the Euler-Mayaruma method to model the discrete nonlinear Schrödinger equation. This mathematical model showcases a phase transition of the order parameter alpha, which was observed numerically. The plots depicted below illustrate the results of the simulation. When the value of $\alpha$ is less than 4, there is no breather present in the system, however, for values of $\alpha$ greater than 4, a stable breather emerges in the system.

<p float="left">
<img src="https://github.com/eurusebr/FeedBack/blob/master/HDNLS(3.7).jpg" alt="alt text" width="400">
<img src="https://github.com/eurusebr/FeedBack/blob/master/HDNLS(4).jpg" alt="alt text" width="400">
<img src="https://github.com/eurusebr/FeedBack/blob/master/HDNLS(5).jpg" alt="alt text" width="400">
<img src="https://github.com/eurusebr/FeedBack/blob/master/HDNLS(8).jpg" alt="alt text" width="400">
</p>

In the plot below you can see the breather amplitude $p_0 = \left | c_n \right |^{2}$ of the breather with $\alpha = 3.7,4,5,8$.
<p align="center">
  <img src="https://github.com/eurusebr/FeedBack/blob/master/p0.jpg" width="450" title="hover text">
</p>




























## References
1. M.Hofmann and B.Drossel, "Dephasing versus collapse: lessons from the tight-binding model with noise", New Journal of Physics, Volume 23, October 2021:
<a href = "https://doi.org/10.1088/1367-2630/ac2ae2"> https://doi.org/10.1088/1367-2630/ac2ae2 </a>

2. M.Hofmann Mastre thesis, "Noise in Tight-BindingSystems with Feedback", July 7,2021
