### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 6cacd9ee-d9fe-4ff3-ba50-c038eebbdb85
begin
	using PlutoUI
	using Plots
	using SparseArrays
    using LinearAlgebra
	using Statistics
	using LaTeXStrings
	using FFTW
end

# ╔═╡ 0e34b19a-5236-48d9-8e39-cfbdac2b484f
using JLD

# ╔═╡ fb2cfca1-4940-4f65-9e8d-429d3aea03e7
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


# ╔═╡ ce30fff2-d65e-11eb-2446-4bd4834c6d15
struct gridData
	dx::Float64
	M::Int64
	
	dt::Float64
	nSteps::Int64
end

# ╔═╡ d683269a-7f87-4ae9-b1a7-47bc9539a30a
grid(p::gridData) = (-p.M:1:p.M ).* p.dx

# ╔═╡ 9e249b64-6806-4a61-ba51-9a46397ad2c4
function initCond(xgrid, sType, μ, σ, p)
	if(sType == 1) #moving gaussian
       return Complex.(1/(σ*sqrt(2*π))^0.5 * exp.(p*im*xgrid - ((xgrid .- μ) .^2 ./ (4*σ^2))))
	end
end

# ╔═╡ 850d7bf5-b3fb-44ba-9c06-352a5bdad489
function GPEsolver(p::gridData, mass::Float64, g::Float64, init, potential)
	
	xgrid = grid(p)
	dx, dt, nSteps = p.dx, p.dt, p.nSteps
	
	n = 2*p.M + 1
	V = potential.(xgrid)
	
	ψ = zeros(Complex, (n, nSteps))
    psi = init
	
	dk = pi/(p.M * p.dx)
	kgrid = dk .* (-p.M:1:p.M)
	
	for j in 1:nSteps
		ψ[:, j] .= psi
		
		psi = psi .* exp.(-0.5 * 1im * dt * (V .+ g * abs2.(psi)))
		psi_k = fftshift(fft(psi)/n)
		psi_k .*= exp.(-0.5 * dt * 1im * (1/mass) * kgrid.^2)
		psi = ifft(ifftshift(psi_k)) * n
		psi = psi .* exp.(-0.5 * 1im * dt * (V .+ g * abs2.(psi)))
		
		psi[1] = psi[end] = 0 + 0im
	end
	
	return (xgrid, ψ)
end

# ╔═╡ cab24696-a0ff-43b9-a6fd-9a7467def699
begin
	p = gridData(0.005, 1600, 0.01, 1000)
	init = initCond(grid(p), 1, 0.0, 1.0, 0.0) .|> Complex
	V(x) = 0.5*x^2 |> Complex
end

# ╔═╡ 1218bc15-9ae7-4a7c-b296-121558605181
md"Solve GPE? $(@bind trigger1 CheckBox())"

# ╔═╡ 5954725a-b81c-4dee-be3b-07a6f3487f71
d = load("./tmp/GroundState.jld")

# ╔═╡ e38ab800-e513-476c-a7b9-ee2f81fb03a8
data = (trigger1) ? GPEsolver(p, 1.0, -5.0, d["psi1"], V) : nothing;

# ╔═╡ c12a7191-e159-4823-8df8-16ca43fdc4ab
function plotProbabilityDensity(xgrid, ψ, m, tskip, xskip, dt)
	#GIF CREATION
    anim = @gif for j = 1:tskip:m
		plot(xgrid[1:xskip:end], abs2.(ψ[1:xskip:end, j]), ylim = (-0.05, 1.0), label = "|\\psi |^2" |> latexstring)
		annotate!([(6.0, 0.35, text("t = $(round(j*dt, digits = 3))", 10))])
    end
end

# ╔═╡ a92047f6-c8d7-4693-9f1c-80e711ee9802
function plotInset(xgrid, ψ, m, tskip, xskip, dt)
	#GIF CREATION
    anim = @gif for j = 1:tskip:m
		plot(xgrid[1:xskip:end], abs2.(ψ[1:xskip:end, j]), ylim = (-0.05, 1.0), label = "|\\psi |^2" |> latexstring)
		annotate!([(6.5, 0.8, text("t = $(round(j*dt, digits = 3))", 10))])
		plot!(xgrid, 0.5*(xgrid).^2, color = :red, inset_subplots = [(1, bbox(0.05, 0, 0.4, 0.4, :top))], subplot=2, framestyle = :box, label = "V(x)", yticks = [], lw = 1.5)
    end
end

# ╔═╡ 46ade0bf-f729-41ce-8473-e0476ee0e0b3
plotInset(data..., p.nSteps, 10, 10, p.dt)

# ╔═╡ c3834557-22b0-41c8-9baf-33b87528d04f
function getGroundState(p::gridData, mass::Float64, g::Float64, init, potential, max_iter)
	
	xgrid = grid(p)
	dx, dt, nSteps = p.dx, p.dt, p.nSteps
	
	n = 2*p.M + 1
	V = potential.(xgrid)
	
    psi = init
	
	prev = psi[100]
	mu_error, mu_old = 1, 1
	count = 0
	
	dk = pi/(p.M * p.dx)
	kgrid = dk .* (-p.M:1:p.M)
	
	while(abs(mu_error) > 1e-8)
		
		psi = psi .* exp.(-0.5 * dt * (V .+ g * abs2.(psi)))
		psi_k = fftshift(fft(psi)/n)
		psi_k .*= exp.(-0.5 * dt * (1/mass) * kgrid.^2)
		psi = ifft(ifftshift(psi_k)) * n
		psi = psi .* exp.(-0.5 * dt * (V .+ g * abs2.(psi)))
		
		psi /= sum(abs2.(psi)) * dx
		
		mu = log(prev/psi[100])/dt
		mu_error = abs(mu - mu_old)/mu
		
		prev = psi[100]
		mu_old = mu
		
		count += 1
		if(count > max_iter) break end
	end
	
	return (count, xgrid, psi)
end

# ╔═╡ 6820554c-0153-4e1c-bd14-129e4b6d71b6
md"Find ground state? $(@bind trigger CheckBox())"

# ╔═╡ 14c6e205-5896-4ac5-86f8-7e7bef0ede27
function plotGroundState(xgrid, psi)
	plot(xgrid, abs2.(psi), ylim = (-0.05, 0.5), label = "")
end

# ╔═╡ 1cc59038-135a-4453-a5df-44fc17ff2fb6
gs = (trigger) ? getGroundState(p, 1.0, 5.0, init, V, 20000) : nothing;

# ╔═╡ 301134a7-06a8-49fb-8d6a-98012d713159
gs[1]

# ╔═╡ 399fa812-4c4b-410d-b222-d34ec7cc0c52
plotGroundState(gs[2:3]...)

# ╔═╡ 4625682b-94d4-4864-bef4-4cf36865d9d6
save("./tmp/GroundState.jld", "psi1", gs[3], "xgrid", gs[2], "g", 0.5, "m", 1, "dt", 1e-3)

# ╔═╡ Cell order:
# ╟─fb2cfca1-4940-4f65-9e8d-429d3aea03e7
# ╠═6cacd9ee-d9fe-4ff3-ba50-c038eebbdb85
# ╠═0e34b19a-5236-48d9-8e39-cfbdac2b484f
# ╠═ce30fff2-d65e-11eb-2446-4bd4834c6d15
# ╠═d683269a-7f87-4ae9-b1a7-47bc9539a30a
# ╠═9e249b64-6806-4a61-ba51-9a46397ad2c4
# ╠═850d7bf5-b3fb-44ba-9c06-352a5bdad489
# ╠═cab24696-a0ff-43b9-a6fd-9a7467def699
# ╟─1218bc15-9ae7-4a7c-b296-121558605181
# ╟─5954725a-b81c-4dee-be3b-07a6f3487f71
# ╠═e38ab800-e513-476c-a7b9-ee2f81fb03a8
# ╟─c12a7191-e159-4823-8df8-16ca43fdc4ab
# ╠═a92047f6-c8d7-4693-9f1c-80e711ee9802
# ╠═46ade0bf-f729-41ce-8473-e0476ee0e0b3
# ╠═c3834557-22b0-41c8-9baf-33b87528d04f
# ╟─6820554c-0153-4e1c-bd14-129e4b6d71b6
# ╠═14c6e205-5896-4ac5-86f8-7e7bef0ede27
# ╠═1cc59038-135a-4453-a5df-44fc17ff2fb6
# ╟─301134a7-06a8-49fb-8d6a-98012d713159
# ╠═399fa812-4c4b-410d-b222-d34ec7cc0c52
# ╠═4625682b-94d4-4864-bef4-4cf36865d9d6
