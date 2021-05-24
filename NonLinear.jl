### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 9442bb6c-adb1-11eb-395c-3d6d225035c4
begin
	using PlutoUI
	using Plots
	using SparseArrays
    using LinearAlgebra
	using Statistics
end

# ╔═╡ 016b11c5-1e58-4e01-979f-1a1c3dbda073
function initCond(xgrid, sType, μ, σ, η)
	
    if(sType == 1)
		return 400 .* η^4 .* exp.(Complex.(0, σ*xgrid)) .* sech.(η*(xgrid .- μ))
	end
	if(sType == 2) #moving gaussian
       return Complex.(1/(σ*sqrt(2*π))^0.5 * exp.(p*im*xgrid - ((xgrid .- μ) .^2 ./ (4*σ^2))))
	end
end

# ╔═╡ 106d3858-3123-487d-9506-c9efd2da0c8c
function schrodingerNonlinear(dx::Float64, xi::Float64, xf::Float64, m::Int64, dt::Float64, k::Float64, wavepacket::Tuple{Float64, Float64, Float64})
    
    xgrid = collect(xi:dx:xf) #discretized spatial grid
    n = length(xgrid)
	
    #INITIAL CONDITIONS
	ψ = zeros(Complex, (n, m))
    ψ[:, 1] = initCond(xgrid, 1, wavepacket...)
	
	#CREATING EVOLUTION MATRIX
	r = dt/(2*dx^2)
	diag = r * ones(Complex, n-1)
	#TIME EVOLUTION LOOP
    for j in 1:(m-1) 
 		α = (k * dt/4) * abs2.(ψ[:, j]) .+ r 
		α₁ = Complex.(1, -α)
		α₂ = Complex.(1, α)
		
		M₁ = LinearAlgebra.SymTridiagonal(α₁, diag)
		M₂ = LinearAlgebra.SymTridiagonal(α₂, -diag)
		
		ψ[:, j+1] = M₁\(M₂*ψ[:, j])
    end
    
    return (xgrid, ψ)
end

# ╔═╡ 0904caee-d44b-49dd-92fb-0656a6fbf4d7
function V(x)
	return (x > -5.0 && x < 5.0) ? 2.5 : 0.0
	#return 0
	#return 0.05*x^2
end

# ╔═╡ 2b51ddd7-7bbc-4110-88d5-5fd524e9287f
function plotSchrodinger(xgrid, ψ, m)
	#GIF CREATION
    anim = @gif for j = 1:10:m
		plot(xgrid, abs2.(ψ[:, j]), ylim = (-0.05, 0.5), label = "")
    end
end

# ╔═╡ 8bb1eb45-3fe6-4dcf-b2a5-930756e8dc4a
m = 2000

# ╔═╡ 52796427-3ce2-4e48-86c1-874c568d18e7
data = schrodingerNonlinear(0.005, -35.0, 35.0, m, 0.01, 1.0, (0.0, 2.0, 0.5));

# ╔═╡ f7bc09b9-2939-4294-8d19-f1f8d5e94fa2
plotSchrodinger(data..., m)

# ╔═╡ Cell order:
# ╠═9442bb6c-adb1-11eb-395c-3d6d225035c4
# ╠═016b11c5-1e58-4e01-979f-1a1c3dbda073
# ╠═106d3858-3123-487d-9506-c9efd2da0c8c
# ╠═0904caee-d44b-49dd-92fb-0656a6fbf4d7
# ╠═2b51ddd7-7bbc-4110-88d5-5fd524e9287f
# ╠═8bb1eb45-3fe6-4dcf-b2a5-930756e8dc4a
# ╠═52796427-3ce2-4e48-86c1-874c568d18e7
# ╠═f7bc09b9-2939-4294-8d19-f1f8d5e94fa2
