### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 8d7029aa-500b-41be-ad01-f40fe5ba8676
begin
	using PlutoUI
	using Plots
	using SparseArrays
    using LinearAlgebra
	using Statistics
	using LaTeXStrings
end

# ╔═╡ 5f44cd92-3e7e-4373-9f8f-bd4ba5f6585c
using FFTW

# ╔═╡ a4a792d8-cebe-11eb-02b6-971141ada5df
function initCond(xgrid, sType, μ, σ, p)
	
    if(sType == 1) #moving gaussian
       return @. Complex(1/(σ*sqrt(2*π))^0.5 * exp(p*im*xgrid - ((xgrid - μ)^2/(4*σ^2))))
	end
end

# ╔═╡ 32b52661-2b7f-4b13-9654-3534ffdc6fb7
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

# ╔═╡ 53fe1015-b1be-422b-be71-8c5a3ef0c111
function V(x)
	return (x > -5.0 && x < -4.5) ? 2.5 : 0.0
	#return 0.0
	#return 0.05*x^2
end

# ╔═╡ 22825407-d967-4f35-9f81-748c805eab66
function plotProbabilityDensity(xgrid, V, ψ, m)
	#GIF CREATION
    anim = @gif for j = 1:10:m
		plot(xgrid, abs2.(ψ[:, j]), ylim = (-0.05, 0.5), label = "")
		plot!(xgrid, V, color = :black, label = "")
    end
end

# ╔═╡ 68d8c951-b72e-4b86-96b2-14712a3e2754
function plotRealImg(xgrid, V, ψ, m)
    anim = @gif for j = 1:10:m
		plot(xgrid, V, color = :black, label = "", ylim = (-0.5, 0.5))
		plot!(xgrid, abs2.(ψ[:, j]), color = :black, label = "|ψ|^2")
		plot!(xgrid, real.(ψ[:, j]), color = :red, label = "Re(ψ)")
		plot!(xgrid, imag.(ψ[:, j]), color = :blue, label = "Im(ψ)")
	end
end

# ╔═╡ c0af67b7-768a-4aff-9e45-f02896b9469e
function plotUncertaintyProduct(dx, xgrid, V, ψ, m)
	n = length(xgrid)
    Mp = SparseArrays.spdiagm(-1 => ones(n-1), 1 => -1*ones(n-1))
	Mp2 = LinearAlgebra.SymTridiagonal(-2*ones(n), 1*ones(n-1))
	xloc = -30.0
	
	anim = @gif for j = 1:10:m
		psi_j = ψ[:, j]
		
		x = dx * sum(xgrid .* abs2.(psi_j)) |> real
		x2 = dx * sum((xgrid.^2) .* abs2.(psi_j)) |> real
	
		p = (-0.5*im) * sum(conj.(psi_j) .* (Mp * psi_j)) |> real
		p2 = (-1/dx) * sum(conj.(psi_j) .* (Mp2 * psi_j)) |> real
		
		sigx = (x2 - x^2) |> sqrt
		sigp = (p2 - p^2) |> sqrt
		
		plot(xgrid, abs2.(psi_j), ylim = (-0.05, 0.5), label = "|ψ|^2")
		plot!(xgrid, V, color = :black, label = "V(x)")
		annotate!([(xloc, 0.45, text("\\sigma_x\\sigma_p = $(round(sigx*sigp, digits = 3))" |> latexstring, 10))])
		annotate!([(xloc, 0.42, text("\\sigma_x = $(round(sigx, digits = 3))" |> latexstring, 10))])
		annotate!([(xloc, 0.39, text("\\sigma_p = $(round(sigp, digits = 3))" |> latexstring, 10))])
		annotate!([(xloc, 0.36, text("\\langle x \\rangle = $(round(x, digits = 3))" |> latexstring, 10))])
		annotate!([(xloc, 0.33, text("\\langle p \\rangle = $(round(-p, digits = 3))" |> latexstring, 10))])
		
    end
end

# ╔═╡ e4ce1d6e-c800-43a9-973d-f92f77032a0b
nSteps = 2000

# ╔═╡ a8376bd0-c839-4367-b943-0ae1c7fcb4e1
data = schrodingerEquation(0.01, -40.0, 30.0, nSteps, 0.025, 5.0, (-25.0, 2.0, 5.5), V);

# ╔═╡ 09bb4b9a-2ba3-4f3a-a1f7-c93bd06f3dc4
plotUncertaintyProduct(0.01, data..., nSteps)

# ╔═╡ 2222f010-a860-4506-be89-1ecffb427231
j = 1500

# ╔═╡ 4e088440-178c-466f-8079-a2b54fa5958e
begin	
	mom_psi = fft(data[3][:, j] .|> ComplexF32) |> fftshift
	pgrid = fftfreq(length(data[1]), data[1][2] - data[1][1] |> (x) -> floor(1/x)) |> fftshift
	
	plot(pgrid, abs2.(mom_psi)/sum(abs2.(mom_psi)))
end

# ╔═╡ Cell order:
# ╠═8d7029aa-500b-41be-ad01-f40fe5ba8676
# ╠═a4a792d8-cebe-11eb-02b6-971141ada5df
# ╠═32b52661-2b7f-4b13-9654-3534ffdc6fb7
# ╠═53fe1015-b1be-422b-be71-8c5a3ef0c111
# ╠═22825407-d967-4f35-9f81-748c805eab66
# ╠═68d8c951-b72e-4b86-96b2-14712a3e2754
# ╠═c0af67b7-768a-4aff-9e45-f02896b9469e
# ╠═e4ce1d6e-c800-43a9-973d-f92f77032a0b
# ╠═a8376bd0-c839-4367-b943-0ae1c7fcb4e1
# ╠═09bb4b9a-2ba3-4f3a-a1f7-c93bd06f3dc4
# ╠═2222f010-a860-4506-be89-1ecffb427231
# ╠═5f44cd92-3e7e-4373-9f8f-bd4ba5f6585c
# ╠═4e088440-178c-466f-8079-a2b54fa5958e
