### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ ebf306ce-d95e-11eb-38a3-33aadcd6a0cb
begin
	using PlutoUI, Plots, LaTeXStrings, FFTW, JLD
	TableOfContents()
end

# ╔═╡ 9fcb961c-7bcb-457d-8bdd-69ee4f34deb2
html"""
<style>

#launch_binder {
	display: none;
}

body.disable_ui main {
		max-width : 95%;
	}

@media screen and (min-width: 1081px) {
	body.disable_ui main {
		margin-left : 10px;
		max-width : 72%;
		align-self: flex-start;
	}
}
</style>
"""

# ╔═╡ e856a3ae-8d28-4d31-975b-0e82ee27eefc
md"""
## Bose-Einstein Condensation

"""

# ╔═╡ 0e04199d-bc24-4c48-882e-a1a70eaefeef
md"""
The phenomenon of Bose-Einstein Condensation occurs for bosonic gases in the low-temperature limit (when the thermal de-broglie wavelength, $\lambda_T$ is of the order of inter-particular distances, $d$ ). 

$$\lambda_T = \frac{\hbar}{\sqrt{2\pi m k_B T}} \approx \left (\frac{V}{N}\right )^{1/3} = d$$

For temperatures below a critical value $T_c$ (or equivalently, above a critical density $n_c$ ), a phase transition occurs in momentum space of the system as a result of which, a macroscopic ($N_0 \approx N$) number of particles settle into the ground state ($k = 0$, momentum state) and behave as a coherent entity (*the condensate*). This exhibits a rich variety of phenomena to study.
"""

# ╔═╡ 9746f8a2-f6df-438d-bdda-da63e61da1fa
md"""
#### Studying the condensate
---
"""

# ╔═╡ 26239009-91ea-4971-894c-74d23b177042
md"""
Ideally, for a system of $N$ particles, we would have to write down the $N$-body hamiltonian and solve the TDSE to obtain the behaviour of the system. However, this is cumbersome, even numerically for $N$ larger than a few hundred particles. There are several approaches to modelling this system, which may involve writing hydrodynamic equations, or a Schrodinger-like time evolution of a quantity.

Here, we describe the condensate with a *"macro-scopic wave function"*, $\Psi$ which also plays the role of the order parameter for this phase transition. We can write the evolution equation for this quantity by using a mean-field approach.
"""

# ╔═╡ 68d93100-2c35-4428-a1cc-f973303ab180
md"""
#### Interactions
---
To begin with, we will consider a simple model of the interactions, under the limit that the bose gas is dilute. We let the particles act as hard balls, and only have contact interaction, i.e; 

$$V(x, x') = g \cdot \delta(x - x')$$

The co-efficient $g$ is dependant on the s-wave scattering co-efficient of the particles, but we will not go into detail here. We simply note that, $g > 0$ creates a repulsive interaction, whereas $g < 0$ creates an attractive interaction.
"""

# ╔═╡ 163544d0-f078-46ed-97f3-ad4f47fced55
md"""
#### The Gross-Pitaevski Equation
---
"""

# ╔═╡ 58f869d5-1bc1-41b2-b99c-4f9033840e71
md"""
The macroscopic wavefunction $\Psi$ will obey an equation similar to the one-body TDSE, with an extra non-linear term arising from the contact interactions. This non-linear term is somewhat representative of the density of the condensate.

$$\left [-\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V(x) + g\cdot|\Psi(x,t)|^2 \right ] \Psi(x,t) = i\hbar\frac{\partial}{\partial t} \Psi(x, t)$$


The quantity $|\Psi(x,t)|^2$ represents the spatial density of the condensate $\equiv n(x, t)$, with the spatial integral giving us the total number of particles in the system. 

$$\int |\Psi(x,t)|^2 dx = \int n(x, t) dx = N$$

We will usually redefine the quantity as $\Psi \equiv \frac{\Psi}{\sqrt{N}}$ so that it normalizes to 1. 
"""

# ╔═╡ d936fe7f-f202-4519-b1ca-c0fda34ebcdb
md"""
### Studying the GPE
---
We first try to find the stationary solutions of this equation. We begin by performing variable separations on the solution; $\Psi(x, \ t) \equiv \psi(x) e^{-i\mu t/\hbar}$, where the quantity $\mu$ is identified as the chemical potential of the system. This gives us the time-independant GPE:


$$\left [-\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V(x) + g\cdot|\psi(x)|^2 \right ] \psi(x) = \mu \psi(x)$$

As this is a non-linear system, the equation is no longer an eigen-value problem. This voids our methods of solving the TISE by diagonalizing the hamiltonian. Instead, we perform another technique that surprisingly works. 

We simulate the **time-dependent GPE in imaginary-time ($t \to -i\tau$)**, and by maintaining normalization, the solution slowly decays into the ground state solution of the time-independent GPE. This can be easily proved for the linear Schrodinger equation, but I have no idea how it works for the non-linear GPE. Regardless, we will proceed with this magical method (and verify that it does indeed provide accurate results).
"""

# ╔═╡ 7987e7b4-39f4-4eca-9ad3-f8008879c822
md"""
## Simulating the GPE
---
A simple finite difference method was attempted but it had serious stability issues, so we proceed instead with a **pseudo-spectral split operator method**. The algorithm proceeds as follows (*stated without any rigorous proofs of the claims made*). 

For the schrodinger equation, we have:

$$i\hbar\frac{\partial}{\partial t} \psi(x, t) = \hat{H}\psi(x,t)$$

where the hamiltonian $\hat{H}$ can be written as $\hat{H} = \hat{T} + \hat{V}$, where $\hat{T} = -\frac{\hbar^2}{2m} \frac{\partial^2}{\partial x^2}$ and $\hat{V} = V(x)$. Since the evolution is unitary and the hamiltonian is time-independent, we can write the evolution in time $\Delta t$ as:

$$\psi(x, t+\Delta t) = e^{-i\Delta t \hat{H}/\hbar}\psi(x,t)$$

Since the operator $T$ and $V$ do not commute, we can approximate the exponent as follows:

$$e^{-i\Delta t \hat{H}/\hbar} \approx e^{-i\Delta t \hat{V}/2\hbar} \cdot e^{-i\Delta t \hat{T}/\hbar} \cdot e^{-i\Delta t \hat{V}/2\hbar} + \mathcal{O}(\Delta t^3)$$

In position space, $V(x)$ is diagonal and simply involves multiplication, and in momentum space, $T(k)$ is diagonal, so we switch between the fourier spaces to apply these operations sequentially like so:

$$\psi(x, t+\Delta t) = e^{-iV(x)\Delta t/2\hbar} \cdot \mathcal{F}^{-1} \left [ e^{-i\hbar k^2\Delta t/2m} \cdot \mathcal{F} \left[ e^{-iV(x)\Delta t/2\hbar} \cdot \psi(x, t) \right ]\right ] + \mathcal{O}(\Delta t^3)$$

Although this method is tailored for the TDSE, it works for the GPE as well, if we replace $V(x) \mapsto V(x) + g|\psi|^2$ and we must utilise the latest value of $\psi$ in the sequential operations. 
"""

# ╔═╡ 542754e5-d4d5-4efb-a58c-b9dc3ee88257
md"""
### Implementation
---
"""

# ╔═╡ f57544b6-26bd-42a3-b60b-74d8fff1e1ae
md"""
We first create a **gridData** structure to hold the specifications of grid spacing and number of iterations. We also define some helper functions to extract data from it.
"""

# ╔═╡ 4d671f92-c422-464c-9068-9260bc91d779
begin
	struct gridData
		dx::Float64
		M::Int64

		dt::Float64
		nSteps::Int64
	end
	
	# Specialized constructor for fixed grid -1.0:dx:1.0
	gridData(dx::Float64, dt::Float64, nSteps::Int64) = gridData(dx, floor(Int, 1/dx),dt, nSteps)
	
	# returns grid points
	grid(p::gridData) = (-p.M:1:p.M).* p.dx
end

# ╔═╡ 618f625a-20c3-4566-951a-f54b1a26ce9d
gridData(1e-3, 1e-2, 1000)

# ╔═╡ 60b23063-e43f-4834-8c4a-1f9a2634cd86
md"""

Now we define a simple initial condition generator, which creates a gaussian wave-packet with some specified momentum.
"""

# ╔═╡ cd0c2f6e-6b6a-4f65-a404-a0fc5e738246
function initCond(xgrid, μ, σ, p)
	
       return Complex.(1/(σ*sqrt(2*π))^0.5 * exp.(p*im*xgrid - ((xgrid .- μ) .^2 ./ (4*σ^2))))
end

# ╔═╡ eda2de55-006c-4fc0-ae53-c36db3ca3c3c
md"""
Here we define the actual solver, by setting up the spatial and fourier space grids as discussed above. 
"""

# ╔═╡ d46ce0c9-cef5-4413-8cf1-3d0379640bd3
function GPEsolver(p::gridData, mass::Float64, g::Float64, init, potential, img_time = false)
	
	# Setting up x-t grid, and gridspecs
	xgrid = grid(p)
	dx, dt, nSteps = p.dx, p.dt, p.nSteps
	
	# Number of spatial grid points
	n = 2*p.M + 1
	V = potential.(xgrid)
	
	ψ = zeros(Complex, (n, nSteps))
    psi = init
	
	# Setting up fourier space grid specs
	dk = pi/(p.M * p.dx)
	kgrid = dk .* (-p.M:1:p.M)
	
	# Running in imaginary time
	if (img_time) dt *= -1im end
		
	# Split-operator time evolution
	for j in 1:nSteps
		ψ[:, j] .= psi
		
		psi = psi .* exp.(-0.5 * 1im * dt * (V .+ g * abs2.(psi)))
		psi_k = fftshift(fft(psi)/n)
		psi_k .*= exp.(-0.5 * dt * 1im * (1/mass) * kgrid.^2)
		psi = ifft(ifftshift(psi_k)) * n
		psi = psi .* exp.(-0.5 * 1im * dt * (V .+ g * abs2.(psi)))
		
		# Re-normalize every time-step
		if (img_time) psi /= (sum(abs2.(psi)) * dx) end
	end
	
	return (xgrid, ψ)
end

# ╔═╡ 799e463d-aeb3-4f53-b3d5-4e134acd5a6b
function GPEsolverFIXEDBC(p::gridData, mass::Float64, g::Float64, init, potential, img_time = false)
	
	# Setting up x-t grid, and gridspecs
	xgrid = grid(p)
	dx, dt, nSteps = p.dx, p.dt, p.nSteps
	
	# Number of spatial grid points
	n = 2*p.M + 1
	V = potential.(xgrid)
	
	ψ = zeros(Complex, (n, nSteps))
    psi = init
	
	# Setting up fourier space grid specs
	dk = pi/(p.M * p.dx)
	kgrid = dk .* (-p.M:1:p.M)
	
	# Running in imaginary time
	if (img_time) dt *= -1im end
		
	# Split-operator time evolution
	for j in 1:nSteps
		ψ[:, j] .= psi
		
		psi = psi .* exp.(-0.5 * 1im * dt * (V .+ g * abs2.(psi)))
		psi_k = fftshift(fft(psi)/n)
		psi_k .*= exp.(-0.5 * dt * 1im * (1/mass) * kgrid.^2)
		psi = ifft(ifftshift(psi_k)) * n
		psi = psi .* exp.(-0.5 * 1im * dt * (V .+ g * abs2.(psi)))
		
		psi[1] = psi[end] = 0 + 0im
		
		# Re-normalize every time-step
		if (img_time) psi /= (sum(abs2.(psi)) * dx) end
	end
	
	return (xgrid, ψ)
end

# ╔═╡ 6523c120-a8a7-460e-860f-6a3e84874bf4
md"""
And some plotting functions to visualize the data.
"""

# ╔═╡ a372aa89-5af4-4c7a-b922-ca462d467020
function plotTimeEvolution(xgrid, ψ, m, tskip, xskip, dt, V)
	#GIF CREATION
    anim = @gif for j = 1:tskip:m
		plot(	xgrid[1:xskip:end], 
				abs2.(ψ[1:xskip:end, j]), 
				ylim = (-0.05, 1.0), 
				label = "|\\psi \\ |^2" |> latexstring)
		
		annotate!([(7.0, 0.8, text("t = $(round(j*dt, digits = 3))", 10))])
		
		plot!(	xgrid, 
				V.(xgrid), 
				color = :red, 
				inset_subplots = [(1, bbox(0.05, 0, 0.4, 0.4, :top))], 
				subplot=2, 
				framestyle = :box, 
				label = "V(x)", 
				lw = 1.5)
    end
end

# ╔═╡ 5a4e1a8a-7840-46be-b2e1-44cd40e78a57
function plotConvergence(xgrid, ψ, m, tskip, xskip, dt)
		mu = Vector{Float64}()
		N = length(xgrid) ÷ 2
	
    anim = @gif for j = 1:tskip:m
		
		p1 = plot(	xgrid[1:xskip:end], 
					abs2.(ψ[1:xskip:end, j]), 
					ylim = (-0.05, 1.0), 
					label = "", 
					title = "Ground state", 
					xlabel = "Position (x)", 
					ylabel = "Condensate density")
		
		annotate!([(3.0, 0.9, text("t = $(round(-j*dt, digits = 2))im" , 10))])
		
		push!(mu, log(abs(ψ[N, j]/ψ[N, j+1])))		
		
		p2 = plot(	1:tskip:j, 
					mu, 
					title = "Convergence", 
					xlabel = "# timesteps", 
					color = :red, 
					lw = 1.5, 
					label = "")
		
		hline!([0], linestyle = :dot, color = :black, label = "")
		
		plot(p1, p2, size = (800, 400))
    end
end

# ╔═╡ a6f5b65c-bad9-4bb7-ae38-2d7208987564
md"""
## Finding the ground state solutions
---
"""

# ╔═╡ fb5d353d-0baa-42d3-9286-006b33012d4f
md"""
##### Particle in a box w/ periodic BC (g > 0): 
"""

# ╔═╡ 3b0ef5f4-4110-4213-8ce8-7424f7330401
let
	p = gridData(0.005, 800, 0.005, 1000)
	init = initCond(grid(p), 0.0, 0.1, 0.0) .|> Complex
	
	V(x) = 0.0 |> Complex
	
	data1 = GPEsolver(p, 1.0, 10.0, init, V, true);
	
	plotConvergence(data1..., 1000, 4, 20, 0.005)
end

# ╔═╡ 4c38dc7f-439f-4b5c-ab80-b3369d140e31
md"""
##### Particle in a box w/ fixed BC (g > 0): 
"""

# ╔═╡ 97b547a5-efb7-40b1-8c5f-91001d535543
let
	p = gridData(0.005, 800, 0.005, 1000)
	init = initCond(grid(p), 0.0, 0.1, 0.0) .|> Complex
	
	V(x) = 0.0 |> Complex
	
	data1 = GPEsolverFIXEDBC(p, 1.0, 10.0, init, V, true);
	
	plotConvergence(data1..., 1000, 4, 5, 0.005)
end

# ╔═╡ 82277f9a-3f3c-4577-bde0-0c41b0f7007e
md"""
##### Particle in a harmonic potential (g > 0): 
"""

# ╔═╡ 20c590e5-65db-4c91-ab06-e129cf8bbc52
let
	p = gridData(0.005, 800, 0.005, 1000)
	init = initCond(grid(p), 0.0, 0.1, 0.0) .|> Complex
	
	V(x) = 0.5 * x^2 |> Complex
	
	data1 = GPEsolver(p, 1.0, 10.0, init, V, true);
	
	plotConvergence(data1..., 1000, 4, 20, 0.005)
end

# ╔═╡ 2843cedd-f63c-4807-b754-84a3972dab78
md"""
### Dynamics in a harmonic potential
---
Let us now take the ground state solution and evolve it in various harmonic potentials to observe collective excitation modes. Let us first solve and store the ground state solution. **We will evolve the state in real-time.**
"""

# ╔═╡ dd053e5e-0fb1-40e7-beb9-157313734899
begin
	# Calculate ground state potential 
	
	params = gridData(0.005, 1600, 0.005, 2000)
	init1 = initCond(grid(params), 0.0, 0.1, 0.0) .|> Complex
	V1(x) = 0.5 * x^2
	
	ψ₀ = GPEsolver(params, 1.0, 10.0, init1, V1, true)[2][:, end];
end

# ╔═╡ 63b76d09-d807-4e0a-803d-97824699869a
md"""
##### Particle in the original harmonic potential (g > 0): 
"""

# ╔═╡ e5cc0452-0170-412f-a121-87e6df44a4f7
let
	V(x) = 0.5 * x^2
	data = GPEsolver(params, 1.0, 10.0, ψ₀, V, false)
	plotTimeEvolution(data..., 2000, 10, 20, 0.005, V)
end

# ╔═╡ fe498e9f-f86b-4578-a94d-847ac23b5335
md"""
##### Particle in the translated harmonic potential (g > 0): 
"""

# ╔═╡ d821f979-4ef8-42eb-bd30-4fb6978e2d4f
let
	V(x) = 0.5 * (x - 2)^2
	data = GPEsolver(params, 1.0, 10.0, ψ₀, V, false)
	plotTimeEvolution(data..., 2000, 10, 20, 0.005, V)
end

# ╔═╡ 5f4edd57-6420-41c9-8e62-3e2b5758227a
md"""
##### Particle in harmonic potential with halved frequency(g > 0): 
"""

# ╔═╡ 8dc92e9f-7672-4e97-a78c-277d64180220
let
	V(x) = 0.25 * x^2
	data = GPEsolver(params, 1.0, 10.0, ψ₀, V, false)
	plotTimeEvolution(data..., 2000, 10, 20, 0.005, V)
end

# ╔═╡ 3db7e1f2-3cbd-40f1-8d0d-f06767595582
md"""
##### Particle in the original harmonic potential (g < 0; attractive): 
"""

# ╔═╡ bd705fad-4000-4142-b0ae-23094ad9adcb
let
	V(x) = 0.5 * x^2
	data = GPEsolver(params, 1.0, -10.0, ψ₀, V, false)
	plotTimeEvolution(data..., 2000, 6, 10, 0.005, V)
end

# ╔═╡ ec56852c-8524-4c0a-a6c3-0ea08619f710
md"""
We note the oscillatory modes in the shifted/varied frequency harmonic traps. In the last simulation, we have flipped the interaction from repulsive to attractive, and we observe a periodic pumping phenomenon where the density peaks in the middle by pulling particles from either side.
"""

# ╔═╡ add6f79c-1c2f-4335-b34d-b59ec3054dad
md"""
#### Implementation note: Units 
---
"""

# ╔═╡ e9ebd213-91d1-4c11-901f-4307d21a9cc9
md"""
Using the actual values of fundamental constants and lengths for the simulation would be unadvisable because these quantities are usually extremely small, which may result in floating point errors. Instead, we set a characteristic length and energy scale for the simulation, and all other quantities are represented in a unit system based off those quantities. 

- For instance, setting  $\hbar = 1$, $m_e = 1$ and $L = 1$ (well width) in the infinite well system uniquely defines a set of units, such that all other quantities are specified as multiples of these. 

- Similarly, in a harmonic trap, the length scale is defined by the trap frequency, $\omega$ like so, $l = \frac{\hbar}{m\omega}$, and the time scale is set as $1/\omega$. 

Without going into the details, these characteristic scales provide a neater simulation setup where we have specified the parameters as simple integers. For this analysis, the specifics not matter much as we are studying only the qualitative behaviour of the time evolution. 
"""

# ╔═╡ 39c73c5e-ddd3-4931-b358-e3b3d8e0ef70
md"""
# References

- *A primer on quantum fluids*; [https://arxiv.org/abs/1605.09580](https://arxiv.org/abs/1605.09580)
"""

# ╔═╡ Cell order:
# ╟─9fcb961c-7bcb-457d-8bdd-69ee4f34deb2
# ╠═ebf306ce-d95e-11eb-38a3-33aadcd6a0cb
# ╟─e856a3ae-8d28-4d31-975b-0e82ee27eefc
# ╟─0e04199d-bc24-4c48-882e-a1a70eaefeef
# ╟─9746f8a2-f6df-438d-bdda-da63e61da1fa
# ╟─26239009-91ea-4971-894c-74d23b177042
# ╟─68d93100-2c35-4428-a1cc-f973303ab180
# ╟─163544d0-f078-46ed-97f3-ad4f47fced55
# ╟─58f869d5-1bc1-41b2-b99c-4f9033840e71
# ╟─d936fe7f-f202-4519-b1ca-c0fda34ebcdb
# ╟─7987e7b4-39f4-4eca-9ad3-f8008879c822
# ╟─542754e5-d4d5-4efb-a58c-b9dc3ee88257
# ╟─f57544b6-26bd-42a3-b60b-74d8fff1e1ae
# ╠═4d671f92-c422-464c-9068-9260bc91d779
# ╠═618f625a-20c3-4566-951a-f54b1a26ce9d
# ╟─60b23063-e43f-4834-8c4a-1f9a2634cd86
# ╠═cd0c2f6e-6b6a-4f65-a404-a0fc5e738246
# ╟─eda2de55-006c-4fc0-ae53-c36db3ca3c3c
# ╠═d46ce0c9-cef5-4413-8cf1-3d0379640bd3
# ╟─799e463d-aeb3-4f53-b3d5-4e134acd5a6b
# ╟─6523c120-a8a7-460e-860f-6a3e84874bf4
# ╠═a372aa89-5af4-4c7a-b922-ca462d467020
# ╠═5a4e1a8a-7840-46be-b2e1-44cd40e78a57
# ╟─a6f5b65c-bad9-4bb7-ae38-2d7208987564
# ╟─fb5d353d-0baa-42d3-9286-006b33012d4f
# ╠═3b0ef5f4-4110-4213-8ce8-7424f7330401
# ╟─4c38dc7f-439f-4b5c-ab80-b3369d140e31
# ╠═97b547a5-efb7-40b1-8c5f-91001d535543
# ╟─82277f9a-3f3c-4577-bde0-0c41b0f7007e
# ╠═20c590e5-65db-4c91-ab06-e129cf8bbc52
# ╟─2843cedd-f63c-4807-b754-84a3972dab78
# ╠═dd053e5e-0fb1-40e7-beb9-157313734899
# ╟─63b76d09-d807-4e0a-803d-97824699869a
# ╠═e5cc0452-0170-412f-a121-87e6df44a4f7
# ╟─fe498e9f-f86b-4578-a94d-847ac23b5335
# ╠═d821f979-4ef8-42eb-bd30-4fb6978e2d4f
# ╟─5f4edd57-6420-41c9-8e62-3e2b5758227a
# ╠═8dc92e9f-7672-4e97-a78c-277d64180220
# ╟─3db7e1f2-3cbd-40f1-8d0d-f06767595582
# ╠═bd705fad-4000-4142-b0ae-23094ad9adcb
# ╟─ec56852c-8524-4c0a-a6c3-0ea08619f710
# ╟─add6f79c-1c2f-4335-b34d-b59ec3054dad
# ╟─e9ebd213-91d1-4c11-901f-4307d21a9cc9
# ╟─39c73c5e-ddd3-4931-b358-e3b3d8e0ef70
