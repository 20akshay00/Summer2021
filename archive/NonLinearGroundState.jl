### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 9442bb6c-adb1-11eb-395c-3d6d225035c4
begin
	using PlutoUI
	using Plots
	using SparseArrays
    using LinearAlgebra
	using Statistics
	using LaTeXStrings
end

# ╔═╡ 4710ef0f-b6e2-4e1c-8115-3ff2307f9e22
using JLD

# ╔═╡ 5ca18933-de71-4c7e-a693-0e467accf88f
md"""
## Gross-Pitaevski equation

"""

# ╔═╡ 016b11c5-1e58-4e01-979f-1a1c3dbda073
function initCond(xgrid, sType, μ, σ, p)
	if(sType == 1) #moving gaussian
       return Complex.(1/(σ*sqrt(2*π))^0.5 * exp.(p*im*xgrid - ((xgrid .- μ) .^2 ./ (4*σ^2))))
	end
end

# ╔═╡ 0904caee-d44b-49dd-92fb-0656a6fbf4d7
function potential(x)
	#return (x > -5.0 && x < 5.0) ? 2.5 : 0.0
	#return 0.0 |> Complex
	return 0.05*x^2 |> Complex
end

# ╔═╡ 106d3858-3123-487d-9506-c9efd2da0c8c
function schrodingerNonlinear(dx::Float64, xi::Float64, xf::Float64, m::Int64, dt::Float64, mass::Float64, g::Float64, wavepacket::Tuple{Float64, Float64, Float64})
    
    xgrid = collect(xi:dx:xf) #discretized spatial grid
    n = length(xgrid)
	V = potential.(xgrid)
	
    #INITIAL CONDITIONS
	ψ = zeros(Complex, (n, m))
    ψ[:, 1] = initCond(xgrid, 1, wavepacket...)
	
	#CREATING EVOLUTION MATRIX
	diag = ones(Complex, n)
	
	α = dt/(2*mass*dx^2) * diag[2: end]
	β = (1 .- dt*(1/(mass*dx^2) .+ V))
	
	#TIME EVOLUTION LOOP
    for j in 1:(m-1) 
		M = LinearAlgebra.SymTridiagonal(β .- dt * g * abs2.(ψ[:, j]), α)
		
		ψ[:, j+1] = M * ψ[:, j]
		
		#normalize
		ψ[:, j+1] /= sqrt(dx * sum(abs2.(ψ[:, j+1])))
    end
    
    return (xgrid, ψ)
end

# ╔═╡ 2b51ddd7-7bbc-4110-88d5-5fd524e9287f
function plotSchrodinger(xgrid, ψ, m, tskip, xskip)
	#GIF CREATION
	mu = Vector{Float64}()
	maxpsi = 1.0
	
    anim = @gif for j = 1:tskip:m
		p1 = plot(xgrid[1:xskip:end], abs2.(ψ[1:xskip:end, j]), ylim = (-0.05, maxpsi), label = "", title = "Ground state", xlabel = "Position (x)", ylabel = "Condensate density")
		annotate!([(3.0, 0.45, text("t = $(-j)im" , 10))])
		
		push!(mu, 100*log(abs(ψ[100, j]/ψ[100, j+1])))		
		p2 = plot(1:tskip:j, mu, title = "Convergence", xlabel = "# timesteps", color = :red, lw = 1.5, label = "")
		hline!([0], linestyle = :dot, color = :black, label = "")
		
		plot(p1, p2, size = (800, 400))
    end
end

# ╔═╡ c570eee9-bdeb-4762-9eba-7af985c969b7
max([1,2,3]...)

# ╔═╡ 8bb1eb45-3fe6-4dcf-b2a5-930756e8dc4a
m = 16000

# ╔═╡ 52796427-3ce2-4e48-86c1-874c568d18e7
data = schrodingerNonlinear(0.005, -4.0, 4.0, m, 0.025, 4000.0, 0.5, (0.0, 1.0, 0.0));

# ╔═╡ f7bc09b9-2939-4294-8d19-f1f8d5e94fa2
plotSchrodinger(data..., m, 100, 20)

# ╔═╡ efc6eab5-b201-4893-a6ae-26aa8805d76d
save("/tmp/GroundState.jld", "psi1", data[2][:, end], "xgrid", data[1])

# ╔═╡ 94d4c67a-d804-4688-b358-5d35a2ce3144
begin	
	mom_psi = fft(data[3][:, j] .|> ComplexF32) |> fftshift
	pgrid = fftfreq(length(data[1]), data[1][2] - data[1][1] |> (x) -> floor(1/x)) |> fftshift
	
	plot(pgrid, abs2.(mom_psi)/sum(abs2.(mom_psi)))
end

# ╔═╡ 85e2b44a-010c-411b-9bbf-654dfe11102e


# ╔═╡ fc4e8efa-8270-41b8-87e6-57565ed26c07
n = 10000

# ╔═╡ 1113575f-bb88-4e65-897f-a9f8985faa2e
function plotSchrodinger2(xgrid, ψ, m, tskip, xskip)
	#GIF CREATION
	mu = Vector{Float64}()
	maxpsi = 0.5
	
    anim = @gif for j = 1:tskip:m
		p1 = plot(xgrid[1:xskip:end], abs2.(ψ[1:xskip:end, j]), ylim = (-0.05, maxpsi), label = "", title = "Ground state", xlabel = "Position (x)", ylabel = "Condensate density")
		annotate!([(3.0, 0.45, text("t = $(j)" , 10))])
    end
end

# ╔═╡ Cell order:
# ╟─9442bb6c-adb1-11eb-395c-3d6d225035c4
# ╟─5ca18933-de71-4c7e-a693-0e467accf88f
# ╠═016b11c5-1e58-4e01-979f-1a1c3dbda073
# ╠═106d3858-3123-487d-9506-c9efd2da0c8c
# ╠═0904caee-d44b-49dd-92fb-0656a6fbf4d7
# ╠═2b51ddd7-7bbc-4110-88d5-5fd524e9287f
# ╠═c570eee9-bdeb-4762-9eba-7af985c969b7
# ╠═8bb1eb45-3fe6-4dcf-b2a5-930756e8dc4a
# ╠═52796427-3ce2-4e48-86c1-874c568d18e7
# ╠═f7bc09b9-2939-4294-8d19-f1f8d5e94fa2
# ╠═4710ef0f-b6e2-4e1c-8115-3ff2307f9e22
# ╠═efc6eab5-b201-4893-a6ae-26aa8805d76d
# ╠═94d4c67a-d804-4688-b358-5d35a2ce3144
# ╠═85e2b44a-010c-411b-9bbf-654dfe11102e
# ╠═fc4e8efa-8270-41b8-87e6-57565ed26c07
# ╠═1113575f-bb88-4e65-897f-a9f8985faa2e
