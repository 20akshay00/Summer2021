### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 7590d87f-ecd6-46e0-9931-a0febbdebc14
begin
	using PlutoUI
	using Plots
	using SparseArrays
    using LinearAlgebra
	using Statistics
end

# ╔═╡ 1ef0df5e-b79f-11eb-3dc2-cf714ea96197
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

# ╔═╡ 81a30c38-4d61-41ae-a654-ebfdf59e7a29
md"## Investigating the relation between TDSE and Heat equation (1D)"

# ╔═╡ 52d10bb6-d20d-4b39-b9e9-d9460db55f0b
md""" 
The TDSE in 1D is given by 

$$i\hbar\frac{\partial \psi(x, \ t)}{\partial t} = \left(-\frac{\hbar^2}{2m} \frac{\partial^2}{\partial x^2} + V(x) \right) \psi(x, \ t)$$ 

If we perform a wick rotation ($$t \to -i\tau$$), the equation now resembles a heat/diffusion equation:

$$\frac{\partial \psi(x, \ i\tau)}{\partial \tau} = \left(\frac{\hbar}{2m} \frac{\partial^2}{\partial x^2} - \frac{1}{\hbar}V(x) \right) \psi(x, \ i\tau)
\hspace{1cm}
\longleftrightarrow
\hspace{1cm}
\frac{\partial \nu(x, \ t)}{\partial t} = \alpha \frac{\partial^2}{\partial x^2} \nu(x, \ t)$$
"""

# ╔═╡ 59691809-7282-46ec-b428-cd8be26d20dd
md"#### Initial Condition"

# ╔═╡ 6a9478e7-1c81-47e5-9266-2d0ec7cbec41
function initCond(xgrid, sType, μ, σ, p)
	
    if(sType == 1) #moving gaussian
       return @. Complex(1/(σ*sqrt(2*π))^0.5 * exp(p*im*xgrid - ((xgrid - μ)^2/(4*σ^2))))
	end
end

# ╔═╡ baf46fd6-511d-4996-acf0-7f4ccf5aef3f
md"#### Crank-Nicholson scheme"

# ╔═╡ 5b25efd6-2392-45ee-b807-37f732b87021
function schrodingerEquationImaginary(dx::Float64, xi::Float64, xf::Float64, m::Int64, dt::Float64, mass::Float64, wavepacket::Tuple{Float64, Float64, Float64}, potential)
    
    xgrid = collect(xi:dx:xf) #discretized spatial grid
    n = length(xgrid)
	
	V = potential.(xgrid)
	
	#CREATING EVOLUTION MATRIX
	diag = ones(Complex, n)
	α = dt/(4*mass*dx^2) * diag[2:end]
	ξ = Complex.(1 .+ (dt/2) * (1/(mass*dx^2) .+ V))
	γ = Complex.(1 .- (dt/2) * (1/(mass*dx^2) .+ V))
	
    M₁ = LinearAlgebra.SymTridiagonal(ξ,-α)
    M₂ = LinearAlgebra.SymTridiagonal(γ, α)
	
    #INITIAL CONDITIONS
	ψ = zeros(Complex, (n, m))
    ψ[:, 1] = initCond(xgrid, 1, wavepacket...)
    
	#TIME EVOLUTION LOOP
    for j in 1:(m-1)
		ψ[:, j+1] = M₁\(M₂*ψ[:, j])
    	
		#normalizing
		ψ[:, j+1] /= sqrt(dx * sum(abs2.(ψ[:, j+1])))
	end
    
    return (xgrid, V, ψ)
end

# ╔═╡ 00656461-1813-4a12-aa1b-53770a8c16a4
md"#### Plotting"

# ╔═╡ 039cd28b-8656-4390-9d08-90a31fb56f5e
function plotProbabilityDensity(xgrid, V, ψ, m)
	#GIF CREATION
    anim = @gif for j = 1:10:m
		plot(xgrid, abs2.(ψ[:, j]), ylim = (-0.05, 0.5), label = "|ψ|^2")
		plot!(xgrid, V, color = :black, label = "V(x)")
    end
end

# ╔═╡ 02787f48-687c-4da6-85c2-f9edb4e02173
md"## Simulation results"

# ╔═╡ a2784ce9-c818-40c3-8625-a7a5575f1513
md"##### Free particle; $$V(x) = 0$$"

# ╔═╡ c0588237-68bc-4982-9479-a4d104115edb
begin
	dataImag1 = schrodingerEquationImaginary(0.01, -30.0, 30.0, 1500, 0.025, 0.5, (0.0, 1.0, 0.0), (x) -> 0.0);
	plotProbabilityDensity(dataImag1..., 1500);
end

# ╔═╡ 6eef0e3a-176d-44a8-83f2-6c381f19e95e
md"This case resembles a diffusion process (with no external drift)."

# ╔═╡ 259912ca-25bc-447c-830c-d41c0c5a5ff3
md"##### Harmonic oscillator; $$V(x) = \frac{1}{2}k x^2$$"

# ╔═╡ 967a285d-db6f-4896-a38b-5d8bac787b66
begin
	dataImag = schrodingerEquationImaginary(0.01, -15.0, 15.0, 1000, 0.025, 0.5, (-10.0, 1.0, 0.0), (x) -> 0.01*x^2);
	plotProbabilityDensity(dataImag..., 1000);
end

# ╔═╡ 698a1b09-93bb-4daa-8d03-9f791e64d37e
md"##### Step potentials inside an infinite well"

# ╔═╡ 020991f3-9f73-4d5f-b063-33ae0037abbd
begin
	dataImag2 = schrodingerEquationImaginary(0.01, -5.0, 5.0, 1000, 0.025, 0.5, (2.0, 1.0, 0.0), (x) -> (x < 0.0) ? 0 : 0.3);
	plotProbabilityDensity(dataImag2..., 1000);
end

# ╔═╡ 8045b71f-f8ef-4a27-8111-94ccc4a7f6e9
begin
	dataImag3 = schrodingerEquationImaginary(0.01, -5.0, 5.0, 1500, 0.025, 0.5, (-2.0, 1.0, 0.0), (x) -> (x > -1.0 && x < 1.0) ? 0.4 : 0.0);
	plotProbabilityDensity(dataImag3..., 1500);
end

# ╔═╡ 44aedc06-3857-42ee-90ee-376e6bde4440
md"#### Observations"

# ╔═╡ 435e68a2-c2fd-4ea3-8d79-55e51d2d0d8c
md"""
- In imaginary time, the solution exponentially decays and no longer preserves normalization. In this implementation,  **we re-scale the wavefunction after time-step** to forcefully preserve the normalization.

"""

# ╔═╡ af841acf-19a8-45e4-8f05-0554837b6c6c
md"
- It is seen that any initial condition of the wavefunction decays to the ground state wavefunction of the given potential."

# ╔═╡ 9f7827a9-5ba0-4638-ba1d-d376ca7c9b25
md"""

###### Calculations
For a bounded system, the hamiltonian H will have a spectrum of bound-state eigen-vectors \{$\psi_k$\} with energy eigen-values \{$E_k$\}. Any initial wave-function can be written as a linear combination of these eigen-vectors which form a complete basis;

$$\phi(x, 0) = \sum_{k}c_k \cdot \psi_k(x)$$

The time evolution can then be written as:

$$\phi(x, t) = \sum_{k}c_k \cdot \psi_k(x) \cdot e^{-iE_kt/\hbar}$$

"""

# ╔═╡ 6457c3b8-94d7-4a8f-a9fc-f666634ebc78
md"""
When we switch to imaginary time ($t \to -i\tau$), the solution goes from oscillatory to exponential decay:

$$\phi(x, \tau) = \sum_{k}c_k \cdot \psi_k(x) \cdot e^{-E_k\tau/\hbar}$$

$$\implies \phi(x, \tau) = e^{-E_0\tau/\hbar}\cdot (c_0\cdot \psi_0(x) + \sum_{k \neq 0}c_k \cdot \psi_k(x) \cdot e^{(E_0-E_k)\tau/\hbar})$$

As a result, $|\psi|^2$ is no longer preserved and we have to re-normalize the wave-function every timestep.

So, asymptotically as $t\to \infty$, any initial condition decays to the ground-state wavefunction of the system.
"""

# ╔═╡ Cell order:
# ╟─1ef0df5e-b79f-11eb-3dc2-cf714ea96197
# ╟─7590d87f-ecd6-46e0-9931-a0febbdebc14
# ╟─81a30c38-4d61-41ae-a654-ebfdf59e7a29
# ╟─52d10bb6-d20d-4b39-b9e9-d9460db55f0b
# ╟─59691809-7282-46ec-b428-cd8be26d20dd
# ╠═6a9478e7-1c81-47e5-9266-2d0ec7cbec41
# ╟─baf46fd6-511d-4996-acf0-7f4ccf5aef3f
# ╠═5b25efd6-2392-45ee-b807-37f732b87021
# ╟─00656461-1813-4a12-aa1b-53770a8c16a4
# ╠═039cd28b-8656-4390-9d08-90a31fb56f5e
# ╟─02787f48-687c-4da6-85c2-f9edb4e02173
# ╟─a2784ce9-c818-40c3-8625-a7a5575f1513
# ╟─c0588237-68bc-4982-9479-a4d104115edb
# ╟─6eef0e3a-176d-44a8-83f2-6c381f19e95e
# ╟─259912ca-25bc-447c-830c-d41c0c5a5ff3
# ╟─967a285d-db6f-4896-a38b-5d8bac787b66
# ╟─698a1b09-93bb-4daa-8d03-9f791e64d37e
# ╟─020991f3-9f73-4d5f-b063-33ae0037abbd
# ╟─8045b71f-f8ef-4a27-8111-94ccc4a7f6e9
# ╟─44aedc06-3857-42ee-90ee-376e6bde4440
# ╟─435e68a2-c2fd-4ea3-8d79-55e51d2d0d8c
# ╟─af841acf-19a8-45e4-8f05-0554837b6c6c
# ╟─9f7827a9-5ba0-4638-ba1d-d376ca7c9b25
# ╟─6457c3b8-94d7-4a8f-a9fc-f666634ebc78
