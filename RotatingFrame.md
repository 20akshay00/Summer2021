$$
\def\bra#1{\left\langle{#1}\right|}
\def\ket#1{\left|{#1}\right\rangle}
\def\braket#1#2{\left\langle{#1}\mid{#2}\right\rangle}
$$

## Schrodinger equation in a rotating frame of reference

The time-dependant Schrodinger equation is given by:
$$
i\hbar\psi_t = \hat{H}\psi
$$
However, based on the frame of reference, the form of the Hamiltonian operator $\hat{H}$ may vary. In an inertial frame we have:
$$
\hat{H} = \frac{\hat{p}^2}{2m} + \hat{V}
$$
A rotating frame of reference is non-inertial, so there may be extra terms appearing in the hamiltonian. 

#### Time evolution in a frame related by a unitary transformation

Let $F \equiv (x, y, z)$ be the lab frame, and $F'\equiv (x',y',z')$ be the rotating frame with constant angular velocity $\vec{\Omega} = \Omega 
\hat{z}$. In the $F$ frame, the evolution of a wave-function $\ket{\psi(t)}$ in a potential is governed by the hamiltonian:
$$
\hat{H} = \frac{\hat{p}^2}{2m} + \hat{V} \hspace{1cm}\rightarrow\hspace{1cm} \hat{H}\ket{\psi(t)} = i\hbar \frac{\partial}{\partial t}\ket{\psi(t)}
$$
If we switch from $F$ to the rotating frame $$F'$$, there exists a unitary transformation $\hat{U(t)}$ such that the wave-function transforms like so, $$\ket{\tilde{\psi}(t)} = \hat{U}(t)\ket{\psi(t)}$$. We can now calculate the new evolution rule in the following way:

$$
i\hbar\frac{\partial}{\partial t}\ket{\tilde{\psi}(t)} = i\hbar\frac{\partial}{\partial t} \left [ \hat{U}(t)\ket{\psi(t)} \right ] = i\hbar \hat{U}(t) \frac{\partial \ket{\psi(t)}}{\partial t} + i\hbar\frac{\partial \hat{U}(t)}{\partial t} \ket{\psi(t)}
$$

$$
\implies i\hbar\frac{\partial}{\partial t}\ket{\tilde{\psi}(t)} = \hat{U}(t)\hat{H}\ket{\psi(t)} + i\hbar \frac{\partial \hat{U}(t)}{\partial t}\ket{\psi(t)}
$$

$$
\implies i\hbar\frac{\partial}{\partial t}\ket{\tilde{\psi}(t)} = \left [\hat{U}(t)\hat{H}\hat{U}^{\dagger}(t) + i\hbar \frac{\partial \hat{U}(t)}{\partial t}\hat{U}^{\dagger}(t) \right ]\ket{\tilde{\psi}(t)}
$$

So, the evolution of the state in frame $F'$ is given by:

$$
i\hbar\frac{\partial}{\partial t}\ket{\tilde{\psi}(t)} = \hat{H'}\ket{\tilde{\psi}(t)} \hspace{1cm}\leftarrow\hspace{1cm} \hat{H'} = \hat{U}\hat{H}\hat{U}^{\dagger} + i\hbar\frac{\partial \hat{U}}{\partial t}\hat{U}^{\dagger}
$$

---

#### Time evolution in a rotating frame
Since angular momentum is the generator of rotations, we can write the unitary operator for the transformation like so:

$$
\hat{U}(t) = \exp\left(\frac{i\Omega t}{\hbar} \hat{L}_z \right)
$$
Where $\hat{L}_z$ is the angular momentum operator along $z$. Since $\hat{L}_z$ commutes with $\hat{H}$, this gives us the new hamiltonian to be:

$$
\hat{H'} = \hat{U}\hat{H}\hat{U}^{\dagger} + i\hbar\frac{\partial \hat{U}}{\partial t}\hat{U}^{\dagger}
$$

$$
\implies \hat{H'} = \hat{H}\hat{U}\hat{U}^{\dagger} + i\hbar\left(\frac{i\Omega \hat{L}_z}{\hbar}\right)\hat{U}\hat{U}^{\dagger}
$$

$$
\implies \hat{H'} = \hat{H} - \Omega\hat{L}_z
$$

More generally, we can write it as: 

$$
\begin{align*}
\hat{H'} &= \hat{H} - \vec{\Omega} \cdot \hat{L} \\
&= \frac{\hat{p}^2}{2m} + \hat{V}(r, \ \theta - \Omega t, \ z) - \vec{\Omega} \cdot \hat{L}
\end{align*}
$$
