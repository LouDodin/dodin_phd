### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ imports
begin
	using PlutoUI
	using Plots
	using Measures
	using Statistics
	using Random
	using Printf
	using BlackBoxOptim
end

# ╔═╡ sliders
md"## Design expérimental — Estimation de γ"

# ╔═╡ slider-ngroups
@bind N_GROUPS Slider(2:1:15, default=5, show_value=true)

# ╔═╡ slider-tmax
@bind T_MAX Slider(10:5:300, default=50, show_value=true)

# ╔═╡ slider-ncells
@bind N_CELLS Slider(10:2:500, default=96, show_value=true)

# ╔═╡ slider-nrep
@bind N_REP Slider(1:1:10, default=3, show_value=true)

# ╔═╡ labels
md"""
| Paramètre | Valeur |
|-----------|--------|
| Nombre de timepoints | $(N_GROUPS) |
| Temps total (jours)  | $(T_MAX) j  |
| Cellules par puit    | $(N_CELLS)  |
| Réplicats            | $(N_REP)    |
| **Cellules totales** | **$(N_GROUPS * N_REP * N_CELLS)** |
"""

# ╔═╡ constants
begin
	const TRUE_r     = 0.574619342477644
	const TRUE_K     = 6.675449070379925e7
	const TRUE_alpha = 2.185e-5
	const TRUE_gamma = 9.47815e-3
	const H0         = 1e5
	const prop_R0    = 1.0 - 1e-15
	const N_MC       = 50
end

# ╔═╡ simulate
function simulate(t_end, gamma, query_times; dt=0.5)
	t_grid = sort(unique(vcat(0.0, collect(0.0:dt:t_end), query_times)))

	S = (1 - prop_R0) * H0
	R = prop_R0 * H0

	pR_grid    = zeros(length(t_grid))
	pR_grid[1] = R / (S + R)

	for i in 2:length(t_grid)
		dt_i  = t_grid[i] - t_grid[i-1]
		H     = S + R
		S_new = max(0.0, S + dt_i * (TRUE_r * S * (1 - H/TRUE_K) - TRUE_alpha * S + gamma * R))
		R_new = max(0.0, R + dt_i * (TRUE_r * R * (1 - H/TRUE_K) + TRUE_alpha * S - gamma * R))
		S, R  = S_new, R_new
		pR_grid[i] = R / (S + R)
	end

	query_idx   = [findfirst(==(t), t_grid) for t in query_times]
	pR_at_query = pR_grid[query_idx]

	return t_grid, pR_grid, pR_at_query
end

# ╔═╡ sample
sample_binomial(n, p) = sum(rand() < p for _ in 1:n) / n

# ╔═╡ fit
function fit_gamma_grid(times, obs)
    γ_vals = exp.(range(log(1e-4), log(0.5), length=500))
    losses = [sum((obs .- propR_at_times(times, γ)).^2) for γ in γ_vals]
    return γ_vals[argmin(losses)]
end

# ╔═╡ monte-carlo
function run_mc(times, true_pR; n_cells, n_rep, n_mc, seed=42)
	Random.seed!(seed)
	gamma_ests   = Float64[]
	gamma_errors = Float64[]

	for _ in 1:n_mc
		obs = [mean(sample_binomial(n_cells, true_pR[i]) for _ in 1:n_rep)
			   for i in eachindex(times)]
		γ_est = fit_gamma(times, obs)
		push!(gamma_ests,   γ_est)
		push!(gamma_errors, abs(γ_est - TRUE_gamma) / TRUE_gamma * 100)
	end

	return gamma_ests, gamma_errors
end

# ╔═╡ compute
begin
	times = [T_MAX * k / N_GROUPS for k in 1:N_GROUPS]
	t_full, pR_full, true_pR = simulate(T_MAX, TRUE_gamma, times)

	gamma_ests, gamma_errors = run_mc(times, true_pR;
		n_cells=N_CELLS, n_rep=N_REP, n_mc=N_MC)

	γ_med   = median(gamma_ests)
	err_med = median(gamma_errors)
	err_p25 = quantile(gamma_errors, 0.25)
	err_p75 = quantile(gamma_errors, 0.75)

	_, pR_fit, _ = simulate(T_MAX, γ_med, times)
end

# ╔═╡ results-table
md"""
## Résultats

| | Valeur |
|---|---|
| γ vrai | $(round(TRUE_gamma, sigdigits=5)) |
| γ estimé (médiane MC) | $(round(γ_med, sigdigits=5)) |
| Erreur médiane | $(round(err_med, digits=1)) % |
| IQR erreur | [$(round(err_p25, digits=1)) % , $(round(err_p75, digits=1)) %] |
"""

# ╔═╡ plot
begin
	eq = TRUE_alpha / (TRUE_alpha + TRUE_gamma)

	p1 = plot(t_full, pR_full,
		lw=4, color=:steelblue,
		label="vraie  γ=$(TRUE_gamma)",
		xlabel="Temps [jours]", ylabel="prop_R [-]",
		title="prop_R(t)",
		ylims=(0.0, 1.05), legend=:bottomleft)

	plot!(p1, t_full, pR_fit,
		lw=2.5, color=:orange, ls=:dash,
		label="fit médian  γ=$(round(γ_med, sigdigits=4))")

	scatter!(p1, times, true_pR,
		ms=8, color=:tomato, markershape=:diamond,
		label="timepoints")

	hline!(p1, [eq], lw=1, ls=:dot, color=:gray,
		label="équilibre ≈ $(round(eq, sigdigits=3))")

	p2 = histogram(gamma_errors,
		bins=15, color=:mediumseagreen, alpha=0.75,
		xlabel="Erreur relative γ (%)", ylabel="Compte",
		title="Distribution erreur γ  (n_mc=$(N_MC))",
		label=false)
	vline!(p2, [err_med], lw=2.5, color=:white,
		label="médiane=$(round(err_med, digits=1))%")
	vline!(p2, [err_p25, err_p75], lw=1.5, color=:lightgray, ls=:dash,
		label="IQR")

	plot(p1, p2,
		layout=(1, 2),
		size=(1400, 500),
		left_margin=10mm, right_margin=6mm,
		top_margin=6mm,   bottom_margin=8mm,
		guidefontsize=13, tickfontsize=11, legendfontsize=10)
end