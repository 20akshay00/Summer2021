---
title: \-
author: Akshay Shankar
date: June 27th
fontsize: 12px
---

# Solving 1D time-independent GPE

$$
\left [ -\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V(x) + g|\psi(x)|^2 \right] \psi(x) = \mu \psi(x)$$



## 

- Run time-dependent GPE in imaginary time ($t \to -i\tau$) till solution converges to ground state.



. . .


$$
\left [ -\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V(x) + g|\psi(x)|^2 \right] \psi(x, t) = i\hbar \frac{\partial}{\partial t} \psi(x, t)
$$

- forward finite difference; stability issues but produced results
- split-operator /pseudo-spectral

---

### Infinite well 

![](/home/akshay/code-repo/Summer2021/ppt/InfiniteWell.gif)

---

### Harmonic potential w/ strong interactions

![](/home/akshay/code-repo/Summer2021/ppt/Harmonic.gif)


---

### Dynamics in harmonic potential

|![Stationary ground state (g > 0)](/home/akshay/code-repo/Summer2021/ppt/Harmonic2.gif){width=85%}| ![Shifted well (g > 0)](/home/akshay/code-repo/Summer2021/ppt/Harmonic1.gif){width=85%} |
| ----------- | ----------- |
| ![Halved well frequency (g > 0)](/home/akshay/code-repo/Summer2021/ppt/Harmonic3.gif){width=85%}     | ![g < 0](/home/akshay/code-repo/Summer2021/ppt/HarmonicAttract.gif){width=85%}     |



# Dipolar BECs

$$U(r) = \frac{C}{4\pi}\frac{1 - 3\cos^2\theta}{r^3}$$

- Anistropic, long range interactions

## 

- Di-polar BECs [https://www.icfo.eu/images/publications/J09-125.pdf](https://www.icfo.eu/images/publications/J09-125.pdf)

- A Primer on Quantum Fluids; [https://arxiv.org/abs/1605.09580](https://arxiv.org/abs/1605.09580)

- Vortex Dynamics, Turbulence and Related Phenomena in Quantum Fluids; [https://www.youtube.com/user/iiptv/](https://youtube.com/playlist?list=PLqTLz9G2bGs2HkIvWusKrJ7kcEp9csnuK)

