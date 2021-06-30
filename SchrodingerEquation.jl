### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ c8576284-9549-11eb-09a9-49c73a096836
begin
	using PlutoUI
	using Plots
	using SparseArrays
    using LinearAlgebra
	using Statistics
end

# ╔═╡ 4dce2ef1-30dd-45db-90f4-c3e8597089a6
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

# ╔═╡ e8a5381b-7810-4239-9b36-8ca4efc68f4e
md"""
# Simulating the 1-D TDSE
"""

# ╔═╡ aad53d7a-eeea-4dd8-bf2d-f012d697dc91
md"""
$$i\hbar\frac{\partial \psi(x, t)}{\partial t} = \left(-\frac{\hbar^2}{2m} \frac{\partial^2}{\partial x^2} + V(x) \right) \psi(x, t)$$ 

"""

# ╔═╡ f10924c1-3f11-4cda-acb2-02aafc605562
md"### Crank-Nicholson scheme"

# ╔═╡ 625ad5be-5f13-4bd5-afbf-d56e4d99b598
md"""
We can write a forward and backward euler expansion for the time derivative (setting $$\hbar = 1$$):

$$\psi(x, t + \Delta t) = \psi(x, t) - i \cdot \hat{H}\psi(x,t) \cdot \Delta t$$

$$\psi(x, t + \Delta t) = \psi(x, t) - i \cdot \hat{H}\psi(x,t + \Delta t) \cdot \Delta t$$

Taking the mean value of these two expressions, we obtain the crank nicholson scheme:

$$\left(1 + \frac{1}{2} i\hat{H}\Delta t\right)\psi(x, t + \Delta t) = \left(1 - \frac{1}{2} i\hat{H}\Delta t\right)\psi(x, t)$$

This gives us a unitary evolution operator, which will preserve the normalization. The second derivative in $$\hat{H}$$ can be approximated using finite differences:

$$\frac{\partial^2 \psi_j^n}{\partial x^2} = \frac{\psi_{j+1}^n - 2\psi_j^n + \psi_{j-1}^n}{(\Delta x)^2}$$

This can then be written in matrix form like so; 

$$M_1 \psi^{n+1} = M_2 \psi^n \hspace{1cm} \implies \hspace{1cm} \psi^{n+1} = M_1 ^{-1} \cdot M_2 \psi^n$$
"""

# ╔═╡ 27f27028-83cb-4b44-86af-435f018169d2
md"""
$$M_1 = \begin{pmatrix}
\xi & -\alpha & 0 & ... & ... & ... & 0\\
-\alpha & \xi & -\alpha & 0 & . &. & .\\
0 & -\alpha & \xi & -\alpha & 0 & . & .\\
\vdots & \ddots & \ddots & \ddots & \ddots & \ddots & \vdots\\
. & . & 0 & -\alpha & \xi & -\alpha & 0\\
. & . & . & 0 & -\alpha & \xi & -\alpha\\
0 & ... & ... & ...& 0 & -\alpha & \xi &
\end{pmatrix}

\hspace{2cm}

M_2 = \begin{pmatrix}
\gamma & \alpha & 0 & ... & ... & ... & 0\\
\alpha & \gamma & \alpha & 0 & . &. & .\\
0 & \alpha & \gamma & \alpha & 0 & . & .\\
\vdots & \ddots & \ddots & \ddots & \ddots & \ddots & \vdots\\
. & . & 0 & \alpha & \gamma & \alpha & 0\\
. & . & . & 0 & \alpha & \gamma & \alpha\\
0 & ... & ... & ...& 0 & \alpha & \gamma &
\end{pmatrix}$$
"""

# ╔═╡ 16911457-303f-4f3a-b64d-5ea72ac3b040
md"""
$$\boxed{\alpha = \frac{i\Delta t}{4m(\Delta x)^2}}
\hspace{2cm}
\boxed{\xi = 1 + \frac{i\Delta t}{2} \cdot \left ( \frac{1}{m(\Delta x)^2} + V(x)\right)}
\hspace{2cm}
\boxed{\gamma = 1 - \frac{i\Delta t}{2} \cdot \left ( \frac{1}{m(\Delta x)^2} + V(x)\right)}$$

"""

# ╔═╡ dafdc30d-2568-490a-8d4b-eec06df17f72
md""" #### Boundary Condition
By virtue of this specific construction of the evolution matrices, the system has fixed boundary conditions ($\psi(x) = 0$ at the boundaries). This effectively emulates an infinite wall on either sides, resulting in a particle-in-a-box system in the absence of any other potential inside.
"""

# ╔═╡ e9d1f155-19a4-40d2-89e7-3b21b7ef856d
md"#### Initial Condition"

# ╔═╡ 208cbf8a-954a-11eb-31c8-014c1b072be0
function initCond(xgrid, sType, μ, σ, p)
	
    if(sType == 1) #moving gaussian
       return @. Complex(1/(σ*sqrt(2*π))^0.5 * exp(p*im*xgrid - ((xgrid - μ)^2/(4*σ^2))))
	end
end

# ╔═╡ fa59ac6f-1e55-4bd7-955e-a71703b54a2b
md"#### Crank-Nicholson evolution"

# ╔═╡ 1fc331e4-954a-11eb-28ba-d17a94789af9
function schrodingerEquation(dx::Float64, xi::Float64, xf::Float64, m::Int64, dt::Float64, mass::Float64, wavepacket::Tuple{Float64, Float64, Float64}, potential)
    
    xgrid = collect(xi:dx:xf) #discretized spatial grid
    n = length(xgrid)
	
	V = potential.(xgrid)
	
	#CONSTRUCTING THE EVOLUTION MATRIX
	diag = ones(Complex, n)
	α = im*dt/(4*mass*dx^2) * diag[2:end]
	ξ = 1 .+ (im * dt/2) * (1/(mass*dx^2) .+ V)
	γ = 1 .- (im * dt/2) * (1/(mass*dx^2) .+ V)
	
    M₁ = LinearAlgebra.SymTridiagonal(ξ,-α)
    M₂ = LinearAlgebra.SymTridiagonal(γ, α)
	
    #INITIAL CONDITIONS
	ψ = zeros(Complex, (n, m))
    ψ[:, 1] = initCond(xgrid, 1, wavepacket...)
    
	#TIME EVOLUTION LOOP
    for j in 1:(m-1)
		ψ[:, j+1] = M₁\(M₂*ψ[:, j])
    end
    
    return (xgrid, V, ψ)
end

# ╔═╡ 9ae10183-b731-4aef-af7a-f048daf4dc92
md"#### Potential function"

# ╔═╡ 270699fa-954a-11eb-3e5b-f1dbdaca0fad
function V(x)
	return (x > -5.0 && x < 5.0) ? 2.5 : 0.0
	#return 0
	#return 0.05*x^2
end

# ╔═╡ af6758f1-302c-401d-9aac-77f8019eea02
md"#### Plotting"

# ╔═╡ 5f14351a-3c56-4d0c-9a79-0d71cbc8eef9
function plotProbabilityDensity(xgrid, V, ψ, m)
	#GIF CREATION
    anim = @gif for j = 1:10:m
		plot(xgrid, abs2.(ψ[:, j]), ylim = (-0.05, 0.5), label = "")
		plot!(xgrid, V, color = :black, label = "")
    end
end

# ╔═╡ d96bc693-7693-4f2e-b119-108bab19cc8b
function plotRealImg(xgrid, V, ψ, m)
    anim = @gif for j = 1:10:m
		plot(xgrid, V, color = :black, label = "", ylim = (-0.5, 0.5))
		plot!(xgrid, abs2.(ψ[:, j]), color = :black, label = "|ψ|^2")
		plot!(xgrid, real.(ψ[:, j]), color = :red, label = "Re(ψ)")
		plot!(xgrid, imag.(ψ[:, j]), color = :blue, label = "Im(ψ)")
	end
end

# ╔═╡ 2aa523ef-178c-4c36-8c1c-a67af7e92ed8
md"## Simulation results"

# ╔═╡ a32b0c6a-895f-4584-90ba-5f17a34a26f2
nSteps = 2200

# ╔═╡ b15fff98-22f9-4689-957f-348dbd5dd02d
data = schrodingerEquation(0.01, -40.0, 40.0, nSteps, 0.025, 5.0, (-25.0, 2.0, 6.0), V);

# ╔═╡ a33fa264-6920-459b-964b-88b453600800
plotProbabilityDensity(data..., nSteps)

# ╔═╡ ca250f60-9b2f-4784-8353-5bc343b1afcd
plotRealImg(data..., nSteps)

# ╔═╡ Cell order:
# ╟─4dce2ef1-30dd-45db-90f4-c3e8597089a6
# ╟─c8576284-9549-11eb-09a9-49c73a096836
# ╟─e8a5381b-7810-4239-9b36-8ca4efc68f4e
# ╟─aad53d7a-eeea-4dd8-bf2d-f012d697dc91
# ╟─f10924c1-3f11-4cda-acb2-02aafc605562
# ╟─625ad5be-5f13-4bd5-afbf-d56e4d99b598
# ╟─27f27028-83cb-4b44-86af-435f018169d2
# ╟─16911457-303f-4f3a-b64d-5ea72ac3b040
# ╟─dafdc30d-2568-490a-8d4b-eec06df17f72
# ╟─e9d1f155-19a4-40d2-89e7-3b21b7ef856d
# ╠═208cbf8a-954a-11eb-31c8-014c1b072be0
# ╟─fa59ac6f-1e55-4bd7-955e-a71703b54a2b
# ╠═1fc331e4-954a-11eb-28ba-d17a94789af9
# ╟─9ae10183-b731-4aef-af7a-f048daf4dc92
# ╠═270699fa-954a-11eb-3e5b-f1dbdaca0fad
# ╟─af6758f1-302c-401d-9aac-77f8019eea02
# ╠═5f14351a-3c56-4d0c-9a79-0d71cbc8eef9
# ╠═d96bc693-7693-4f2e-b119-108bab19cc8b
# ╟─2aa523ef-178c-4c36-8c1c-a67af7e92ed8
# ╠═a32b0c6a-895f-4584-90ba-5f17a34a26f2
# ╠═b15fff98-22f9-4689-957f-348dbd5dd02d
# ╠═a33fa264-6920-459b-964b-88b453600800
# ╠═ca250f60-9b2f-4784-8353-5bc343b1afcd
