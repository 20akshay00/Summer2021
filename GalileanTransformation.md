## Galilean transformation of the Schrodinger Equation

We would expect that the form of the TDSE is preserved under a galilean  transformation. Let us see whether we can arrive at this result starting from the equation in one frame, and performing the transformations.

---

We'll work in 1D for simplicity. Let us take a frame of reference $F$ with co-ordinates $x$ and $t$, and another frame $F'$ with co-ordinates $x'$ and $t'$ such that $F'$ is moving at a constant velocity $v$ relative to $F$. The co-ordinates are related by a galilean boost like so:
$$
x' = x + vt  \hspace{2cm} t' = t
$$
The potential energy is given by $V(x, t)$ in $F$, and by $V'(x,t)$ in $F'$, which are related in the following way:

$$
V'(x,t) = V(x', t') \equiv V(x + vt, t)
$$
In $F$, the TDSE has the form;
$$
\boxed{\left [ \frac{-\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V(x, t) \right ] \psi(x, t)= i\hbar\frac{\partial }{\partial t} \psi(x, t)}
$$
where, $\psi(x, t)$ is the wavefunction in $F$.  On the other hand, in frame $F'$, the wavefunction may have a different form $\psi'(x, t)$ and if the TDSE is galilean covariant, we would expect it to satisfy the equation:
$$
\boxed{\left [ \frac{-\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V'(x, t) \right ] \psi'(x, t)= i\hbar\frac{\partial }{\partial t} \psi'(x, t)}
$$
Since the jacobian of this transformation is 1, we know that the probability density at a specific point in space must be the same in $F$ and $F'$; 
$$
|\psi'(x, t)|^2 = |\psi(x',t')|^2
$$

$$
\implies \boxed{\psi'(x,t) = e^{-if(x',t')}\psi(x',t')}
$$

We can now find the transformation of the derivatives under the this galilean boost:
$$
\frac{\partial}{\partial x} = \frac{\partial x'}{\partial x}\frac{\partial}{\partial x'} + \frac{\partial t'}{\partial x}\frac{\partial}{\partial t'} = \frac{\partial}{\partial x'}
\\
\frac{\partial}{\partial t} = \frac{\partial x'}{\partial t}\frac{\partial}{\partial x'} + \frac{\partial t'}{\partial t}\frac{\partial}{\partial t'} = \frac{\partial}{\partial t'} + v \frac{\partial}{\partial x'}
$$

We still do not know the form of f(x', t'),  but if the TDSE is galilean covariant, the consistency of the equations in the two frames will uniquely determine a form of this function, proving that the equation is indeed covariant under this transformation.

### Finding the form of  $\mathbf{f(x, t)}$

We begin with the TDSE in the frame $F'$:
$$
\left [ \frac{-\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V'(x, t) \right ] \psi'(x, t)= i\hbar\frac{\partial }{\partial t} \psi'(x, t)
$$
Substituting the relations between the new derivatives and wavefunctions:
$$
\left [ \frac{-\hbar^2}{2m}\frac{\partial^2}{\partial x'^2} + V(x', t') \right ] e^{-if(x', t')}\psi(x', t')= i\hbar\left[\frac{\partial }{\partial t'} + v \frac{\partial}{\partial x'}\right ]e^{-if(x', t')}\psi(x', t')
$$
We will now drop the labels (prime) for ease of writing, since all quantities are represented in terms of $(x', t')$ which can just as well be written as the variables $(x, t)$.

$$
\equiv \left [ \frac{-\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V(x, t) \right ] e^{-if(x,t)}\psi(x, t)= i\hbar\left[\frac{\partial }{\partial t} + v\frac{\partial }{\partial x} \right] e^{-if(x,t)}\psi(x, t)
$$

$$
\frac{-\hbar^2}{2m} \left [ (\psi_x - if'\psi)\cdot e^{-if} \right ]' + V\psi \cdot e^{-if} = i\hbar\left[\psi_t - if_t\psi + v\psi_x - ivf_x\psi \right] \cdot e^{-if}
$$

$$
\frac{-\hbar^2}{2m} \left [\psi_{xx} - 2if_x\psi_x + (f_x^2 - if_{xx})\psi \right] \cdot \cancel{e^{-if}} + V\psi \cdot \cancel{e^{-if}} = i\hbar \left[\psi_t - if_t\psi + v\psi_x - ivf_x\psi \right] \cdot \cancel{e^{-if}}
$$

$$
\frac{-\hbar^2}{2m}\psi_{xx} + V\psi + \left[-\frac{\hbar^2}{2m}(-2if_x\psi_x + (f_x^2-if_{xx})\psi)\right] = i\hbar\psi_t + i\hbar\left[ -if_t\psi + v\psi_x - ivf_x\psi\right]
$$

$$
\boxed{\underbrace{\frac{-\hbar^2}{2m}\psi_{xx} + V\psi = i\hbar\psi_t}_{\text{TDSE in frame F}} - \underbrace{\left [ -\frac{\hbar^2}{2m}f_x^2 + \frac{i\hbar}{2m}f_{xx} - \hbar f_t - \hbar vf_x\right]\psi - \left[\frac{i\hbar^2}{m}f_x - i\hbar v \right]\psi_x}_{\text{must vanish for arbitrary }\psi}}
$$

---

To preserve the form of the TDSE in frame F, the last two terms must vanish for arbitrary $\psi$.

$$
\implies -\frac{\hbar^2}{2m}f_x^2 + \frac{i\hbar}{2m}f_{xx} - \hbar f_t - \hbar vf_x = 0 \hspace{2cm} \frac{i\hbar^2}{m}f_x - i\hbar v = 0
$$


The second equation gives us that $f(x,t) = mvx/\hbar + g(t)$. This sets $f_{xx} = 0$. Substituting this in equation (1) gives us:

$$
\frac{-\hbar^2}{2m}\cdot \frac{m^2v^2}{\hbar^2} + 0 - \hbar f_t -mv^2 = 0
$$

$$
f_t = \frac{mv^2}{2\hbar} \implies f(x,t) = g(x) + \frac{mv^2 t}{\hbar}
$$

This gives us the expression for $f$:
$$
f(x,t) = \frac{1}{\hbar}\left[mvx - \frac{1}{2}mv^2t\right]
$$


---

So we see that if $\psi'(x,t)$ in frame $F'$ satisfies the equation:

$$
\left [ \frac{-\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V'(x, t) \right ] \psi'(x, t)= i\hbar\frac{\partial }{\partial t} \psi'(x, t)
$$


Then $\psi(x,t)$ in frame $F$ satisfies the equation: 

$$
\left [ \frac{-\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V(x, t) \right ] \psi(x, t)= i\hbar\frac{\partial }{\partial t} \psi(x, t)
$$
where the $\psi\to\psi'$ transforms like so:
$$
\psi'(x,t) = e^{-\frac{i}{\hbar}\left[mv(x+vt) - \frac{1}{2}mv^2t\right]}\cdot \psi(x+vt,t)
$$


The existence of a solution $f(x,t)$ for imposing this form invariance of the TDSE indicates the galilean covariance of this equation.

---

