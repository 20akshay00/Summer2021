### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ b437fb1a-f1a0-4be8-949d-23bfdd4827ce
html"<style>
main {
    max-width:75%;
    padding: unset;
    margin-right: unset !important;
    margin-bottom: 20px;
    align-self: center !important;
}
pluto-helpbox {
    visibility: hidden;
}
</style>"

# ╔═╡ 4a1f77e0-bc80-11eb-3c3b-931ce8024e18
md"### Galilean transformation of the Schrodinger Equation"

# ╔═╡ 4b609a93-607e-4ab5-89c9-17225412d10c
md"""
We would expect that the form of the TDSE is invariant under a galilean transformation. Let us see whether we can arrive at this result starting from the equation in one frame, and performing the transformations.
"""

# ╔═╡ 0b0a4161-15f6-41f5-b4d1-ecb43b08b828
md"""

---


"""

# ╔═╡ b32b5855-b7ab-4eef-8ae9-78d958f089c8
md"""
We'll work in 1D for simplicity. Let us take a frame of reference $F$ with co-ordinates $x$ and $t$, and another frame $F'$ with co-ordinates $x'$ and $t'$ such that $F'$ is moving at a constant velocity $v$ relative to $F$. The co-ordinates are related by a galilean boost like so:

$$x = x' + vt'  \hspace{2cm} t = t'$$

The potential energy is given by $V(x, t)$ in $F$, and by $V'(x',t')$ in $F'$, with the assumption;

$$V'(x',t') = V(x, t)$$

In $F'$, the TDSE has the form 

$$\left [ \frac{-\hbar^2}{2m}\frac{\partial^2}{\partial x'^2} + V'(x', t') \right ] \psi'(x', t')= i\hbar\frac{\partial }{\partial t'} \psi'(x', t')$$

where, $\psi'(x', t')$ is the wavefunction in $F'$. In $F$, the wavefunction is $\psi(x, t)$ and we find the equation it satisfies by transforming the operators in the above equation. 

$$\frac{\partial}{\partial x'} = \frac{\partial x}{\partial x'}\frac{\partial}{\partial x} + \frac{\partial t}{\partial x'}\frac{\partial}{\partial t} = \frac{\partial}{\partial x}$$

$$\frac{\partial}{\partial t'} = \frac{\partial x}{\partial t'}\frac{\partial}{\partial x} + \frac{\partial t}{\partial t'}\frac{\partial}{\partial t} = \frac{\partial}{\partial t} + v \frac{\partial}{\partial x}$$

We also have another piece of information; since the jacobian of this transformation is 1, we know that the probability density at a specific point must be the same in $F$ and $F'$; 

$$|\psi(x, t)|^2 = |\psi'(x',t')|^2$$

$$\implies \psi(x,t) = e^{if(x,t)}\psi'(x',t')$$

With this, we can substitute these new operators and wavefunctions into the TDSE in $F'$.
"""

# ╔═╡ 3dc47f24-51a7-496b-8be0-44b00ed1d199
md"""

---

"""

# ╔═╡ 9fdecb54-979b-42ed-84c5-aa64409b0cae
md"##### Substituting and re-arranging terms"

# ╔═╡ 649f2a36-db54-400d-a38b-3986383ea681
md"""

$$\left [ \frac{-\hbar^2}{2m}\frac{\partial^2}{\partial x'^2} + V'(x', t') \right ] \psi'(x', t')= i\hbar\frac{\partial }{\partial t'} \psi'(x', t')$$


$$\left [ \frac{-\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V(x, t) \right ] e^{-if(x,t)}\psi(x, t)= i\hbar\left[\frac{\partial }{\partial t} + v\frac{\partial }{\partial x} \right] e^{-if(x,t)}\psi(x, t)$$

$$\frac{-\hbar^2}{2m} \left [ (\psi_x - if'\psi)e^{-if} \right ]' = i\hbar e^{-if}\left[\psi_t - if_t\psi + v\psi_x - ivf_x\psi \right]$$

$$\frac{-\hbar^2}{2m} \left [\psi_{xx} - 2if_x\psi_x + (f_x^2 - if_{xx})\psi \right]e^{-if}= i\hbar e^{-if}\left[\psi_t - if_t\psi + v\psi_x - ivf_x\psi \right]$$

$$\frac{-\hbar^2}{2m}\psi_{xx} + V\psi + \left[-\frac{\hbar^2}{2m}(-2if_x\psi_x + (f_x^2-if_{xx})\psi)\right] = i\hbar\psi_t + i\hbar\left[ -if_t\psi + v\psi_x - ivf_x\psi\right]$$

$$\frac{-\hbar^2}{2m}\psi_{xx} + V\psi = i\hbar\psi_t - \left [ -\frac{\hbar^2}{2m}f_x^2 + \frac{i\hbar}{2m}f_{xx} - \hbar f_t - \hbar vf_x\right]\psi - \left[\frac{i\hbar^2}{m}f_x - i\hbar v \right]\psi_x$$

---

"""

# ╔═╡ 10ca7cd9-47f6-4e1f-ae01-1455622ac02b
md"""
To obtain the form of the TDSE that we want, the last two terms must go to 0 for arbitrary $\psi$.

$$\implies -\frac{\hbar^2}{2m}f_x^2 + \frac{i\hbar}{2m}f_{xx} - \hbar f_t - \hbar vf_x = 0 \hspace{2cm} \frac{i\hbar^2}{m}f_x - i\hbar v = 0$$

The second equation gives us that $f(x,t) = mvx/\hbar + g(t)$. This sets $f_{xx} = 0$. Substituting this in equation (1) gives us:

$$\frac{-\hbar^2}{2m}\cdot \frac{m^2v^2}{\hbar^2} + 0 - \hbar f_t -mv^2 = 0$$

$$f_t = \frac{mv^2}{2\hbar} \implies f(x,t) = g(x) + \frac{mv^2 t}{\hbar}$$

This gives us the expression for $f$:

$$f(x,t) = \frac{1}{\hbar}\left[mvx - \frac{1}{2}mv^2t\right]$$

---

"""

# ╔═╡ cdb9d6b5-54fc-4cd2-9063-ae8b04ad3eff
md"""
So we see that if $\psi'(x',t')$ in frame $F'$ satisfies the equation:

$$\left [ \frac{-\hbar^2}{2m}\frac{\partial^2}{\partial x'^2} + V'(x', t') \right ] \psi'(x', t')= i\hbar\frac{\partial }{\partial t'} \psi'(x', t')$$

Then $\psi(x,t)$ in frame $F$ satisfies the equation: 

$$\left [ \frac{-\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V(x, t) \right ] \psi(x, t)= i\hbar\frac{\partial }{\partial t} \psi(x, t)$$

where the $\psi\to\psi'$ transforms like so:

$$\psi'(x',t') = e^{-\frac{i}{\hbar}\left[mvx - \frac{1}{2}mv^2t\right]}\cdot \psi(x,t)$$

The existence of a solution $f(x,t)$ for imposing this form invariance of the TDSE indicates the galilean invariance of this equation.
"""

# ╔═╡ Cell order:
# ╟─b437fb1a-f1a0-4be8-949d-23bfdd4827ce
# ╟─4a1f77e0-bc80-11eb-3c3b-931ce8024e18
# ╟─4b609a93-607e-4ab5-89c9-17225412d10c
# ╟─0b0a4161-15f6-41f5-b4d1-ecb43b08b828
# ╟─b32b5855-b7ab-4eef-8ae9-78d958f089c8
# ╟─3dc47f24-51a7-496b-8be0-44b00ed1d199
# ╟─9fdecb54-979b-42ed-84c5-aa64409b0cae
# ╟─649f2a36-db54-400d-a38b-3986383ea681
# ╟─10ca7cd9-47f6-4e1f-ae01-1455622ac02b
# ╟─cdb9d6b5-54fc-4cd2-9063-ae8b04ad3eff
