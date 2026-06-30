using Plots
using Measures
using LaTeXStrings
using Statistics
using Random
using Printf
using BlackBoxOptim
using Distributions

# ─────────────────────────────────────────────────────────────────────────────
# PARAMÈTRES VRAIS (simulation)
# ─────────────────────────────────────────────────────────────────────────────
const TRUE_r     = 0.574619342477644
const TRUE_K     = 6.675449070379925e7
const TRUE_alpha = 2.185e-5
const TRUE_gamma = 9.47815e-3
const H0         = 1e5
const prop_R0    = 1.0 - 1e-15

# ─────────────────────────────────────────────────────────────────────────────
# SIMULATION ODE
# ─────────────────────────────────────────────────────────────────────────────
function simulate(t_end, gamma, query_times; dt=0.5)

    base_grid = collect(0.0:dt:t_end)
    t_grid = sort(unique(vcat(base_grid, query_times)))

    time_index = Dict(t => i for (i, t) in enumerate(t_grid))

    S = (1 - prop_R0) * H0
    R = prop_R0 * H0

    pR_grid = zeros(length(t_grid))
    pR_grid[1] = R / (S + R)

    for i in 2:length(t_grid)
        dt_i = t_grid[i] - t_grid[i-1]
        H = S + R

        S_new = max(0.0, S + dt_i * (
            TRUE_r * S * (1 - H/TRUE_K) -
            TRUE_alpha * S +
            gamma * R
        ))

        R_new = max(0.0, R + dt_i * (
            TRUE_r * R * (1 - H/TRUE_K) +
            TRUE_alpha * S -
            gamma * R
        ))

        S, R = S_new, R_new
        pR_grid[i] = R / (S + R)
    end

    idx = [time_index[t] for t in query_times]

    return t_grid, pR_grid, isempty(idx) ? Float64[] : pR_grid[idx]
end

# ─────────────────────────────────────────────────────────────────────────────
# ÉCHANTILLONNAGE BINOMIAL STABLE
# ─────────────────────────────────────────────────────────────────────────────
sample_binomial(n, p, rng) = rand(rng, Binomial(n, p)) / n

# ─────────────────────────────────────────────────────────────────────────────
# FIT γ
# ─────────────────────────────────────────────────────────────────────────────
function fit_gamma(times, obs)

    function objective(x)
        _, _, pred = simulate(maximum(times), x[1], times)
        return sum((obs .- pred).^2)
    end

    res = bboptimize(
        objective;
        SearchRange    = [(1e-5, 0.5)],
        NumDimensions  = 1,
        Method         = :adaptive_de_rand_1_bin_radiuslimited,
        PopulationSize = 200,
        MaxSteps       = 80_000,
        TraceMode      = :silent
    )

    return best_candidate(res)[1]
end

# ─────────────────────────────────────────────────────────────────────────────
# PARAMÈTRES EXPÉRIMENTAUX
# ─────────────────────────────────────────────────────────────────────────────
const T_MAX   = 20
const N_CELLS = 96
const N_REP   = 5
const N_MC    = 300
const N_TIME  = 1

Random.seed!(42)

t_full, pR_full, _ = simulate(T_MAX, TRUE_gamma, Float64[])
times = [T_MAX/N_TIME*i for i in 1:N_TIME]

_, _, true_pR = simulate(T_MAX, TRUE_gamma, times)

# ─────────────────────────────────────────────────────────────────────────────
# MONTE CARLO
# ─────────────────────────────────────────────────────────────────────────────
γ_est   = zeros(N_MC)
γ_error = zeros(N_MC)
all_obs = Vector{Vector{Float64}}(undef, N_MC)

rngs = [MersenneTwister(1234 + i) for i in 1:N_MC]

Threads.@threads for i in 1:N_MC
    rng = rngs[i]

    obs = [
        mean(sample_binomial(N_CELLS, true_pR[j], rng) for _ in 1:N_REP)
        for j in eachindex(times)
    ]

    all_obs[i] = obs

    γ_est[i] = fit_gamma(times, obs)
    γ_error[i] = abs(γ_est[i] - TRUE_gamma) / TRUE_gamma * 100

    @printf("MC %3d/%d γ=%.5f err=%.2f%%\n", i, N_MC, γ_est[i], γ_error[i])
end

# ─────────────────────────────────────────────────────────────────────────────
# STATISTIQUES
# ─────────────────────────────────────────────────────────────────────────────
γ_med = median(γ_est)

err_med = median(γ_error)
err_p25 = quantile(γ_error, 0.25)
err_p75 = quantile(γ_error, 0.75)

println("\n=== Résultats ===")
println("γ vrai   = ", TRUE_gamma)
println("γ médian = ", γ_med)
println("erreur médiane = ", err_med)

# ─────────────────────────────────────────────────────────────────────────────
# COURBE FIT MÉDIANE
# ─────────────────────────────────────────────────────────────────────────────
_, _, _ = simulate(T_MAX, γ_med, times)

# ─────────────────────────────────────────────────────────────────────────────
# PLOT
# ─────────────────────────────────────────────────────────────────────────────
p1 = plot(t_full, pR_full,
          lw=4, color=:steelblue,
          label="vraie γ=$(TRUE_gamma)",
          xlabel="Temps [jours]", ylabel="prop_R [-]",
          title="prop_R(t) — T=$(Int(T_MAX)) j | N=$(N_CELLS) cellules | $(N_TIME) timepoints | $(N_REP) réplicats",
          ylims=(0.0, 1.05), legend=:bottomleft,
          titlefontsize=14)

for obs in all_obs
    scatter!(p1, times, obs, ms=4, color=:gray, alpha=0.15, label=false)
end

scatter!(p1, times, true_pR,
         ms=8, color=:tomato, markershape=:diamond,
         label="vraie aux timepoints")

p2 = histogram(γ_error,
               bins=15, color=:mediumseagreen, alpha=0.75,
               xlabel="Erreur relative γ (%)", ylabel="Compte",
               title="Distribution erreur γ (n_mc=$(N_MC))",
               label=false)

vline!(p2, [err_med], lw=2.5, color=:white,
       label="médiane=$(round(err_med, digits=1))%")

vline!(p2, [err_p25, err_p75], lw=1.5, color=:lightgray, ls=:dash,
       label="IQR")

pl = plot(p1, p2,
          layout=(1, 2),
          size=(3000, 1200),
          left_margin=15mm, right_margin=10mm,
          top_margin=10mm, bottom_margin=12mm,
          plot_title="Design : T=$(Int(T_MAX)) j | N=$(N_CELLS) | $(N_TIME) timepoints | $(N_REP) réplicats",
          plot_titlefontsize=20,
          guidefontsize=16, tickfontsize=14, legendfontsize=13)

savefig(pl, joinpath(@__DIR__, "gamma_single_design.png"))
display(pl)

println("\nFigure sauvegardée : gamma_single_design.png")