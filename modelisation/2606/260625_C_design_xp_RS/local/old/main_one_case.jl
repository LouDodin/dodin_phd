using Plots
using Measures
using LaTeXStrings
using Statistics
using Random
using Printf
using BlackBoxOptim

# ── Paramètres vrais ─────────────────────────────────────────────────────────
const TRUE_r     = 0.574619342477644
const TRUE_K     = 6.675449070379925e7
const TRUE_alpha = 2.185e-5
const TRUE_gamma = 9.47815e-3
const H0         = 1e5
const prop_R0    = 1.0 - 1e-15

# ── Simulation ODE ────────────────────────────────────────────────────────────
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

# ── Échantillonnage binomial (thread-safe RNG) ───────────────────────────────
sample_binomial(n, p, rng) = sum(rand(rng) < p for _ in 1:n) / n

# ── Fitting γ (BlackBoxOptim) ─────────────────────────────────────────────────
function fit_gamma(times, obs)
    function objective(x)
        _, _, pR_pred = simulate(maximum(times), x[1], times)
        return sum((obs .- pR_pred).^2)
    end

    res = bboptimize(
        objective;
        SearchRange          = [(1e-4, 0.5)],
        NumDimensions        = 1,
        Method               = :adaptive_de_rand_1_bin_radiuslimited,
        PopulationSize       = 250,
        MaxSteps             = 200_000,
        DifferentialWeight   = 0.9,
        CrossoverProbability = 0.9,
        TraceMode            = :silent,
    )
    return best_candidate(res)[1]
end

# ── Paramètres expérimentaux ──────────────────────────────────────────────────
const N_GROUPS = 5
const T_MAX    = 35.0
const N_CELLS  = 32
const N_REP    = 3
const N_MC     = 500

Random.seed!(42)

times = [T_MAX * k / N_GROUPS for k in 1:N_GROUPS]
t_full, pR_full, true_pR = simulate(T_MAX, TRUE_gamma, times)

@printf("Timepoints : %s jours\n", join(round.(times, digits=1), ", "))
@printf("prop_R vrai aux timepoints : %s\n\n",
        join(round.(true_pR, sigdigits=3), ", "))

# ── Préparation MC ────────────────────────────────────────────────────────────
gamma_ests   = Vector{Float64}(undef, N_MC)
gamma_errors = Vector{Float64}(undef, N_MC)
all_obs      = Vector{Vector{Float64}}(undef, N_MC)

rngs = [MersenneTwister(42 + i) for i in 1:N_MC]

# ── Monte Carlo parallélisé ──────────────────────────────────────────────────
Threads.@threads for mc in 1:N_MC
    rng = rngs[mc]

    obs = [
        mean(sample_binomial(N_CELLS, true_pR[i], rng) for _ in 1:N_REP)
        for i in eachindex(times)
    ]

    all_obs[mc] = obs

    γ_est = fit_gamma(times, obs)
    gamma_ests[mc] = γ_est
    gamma_errors[mc] = abs(γ_est - TRUE_gamma) / TRUE_gamma * 100

    @printf("MC %2d/%d — γ_est=%.5f  erreur=%.1f%%\n",
            mc, N_MC, γ_est, gamma_errors[mc])
end

γ_med   = median(gamma_ests)
err_med = median(gamma_errors)
err_p25 = quantile(gamma_errors, 0.25)
err_p75 = quantile(gamma_errors, 0.75)

@printf("\n=== Résultats ===\n")
@printf("γ vrai          = %.6f\n", TRUE_gamma)
@printf("γ estimé (med)  = %.6f\n", γ_med)
@printf("Erreur médiane  = %.1f%%\n", err_med)
@printf("IQR erreur      = [%.1f%%, %.1f%%]\n", err_p25, err_p75)

# ── Courbe fittée médiane ─────────────────────────────────────────────────────
_, pR_fit, _ = simulate(T_MAX, γ_med, times)

# ── Plot ──────────────────────────────────────────────────────────────────────
p1 = plot(t_full, pR_full,
          lw=4, color=:steelblue,
          label="vraie  γ=$(TRUE_gamma)",
          xlabel="Temps [jours]", ylabel="prop_R [-]",
          title="prop_R(t)  —  T=$(Int(T_MAX)) j | N=$(N_CELLS) cellules | $(N_GROUPS) timepoints | $(N_REP) réplicats",
          ylims=(0.0, 1.05), legend=:bottomleft,
          titlefontsize=14)

for obs in all_obs
    scatter!(p1, times, obs, ms=4, color=:gray, alpha=0.15, label=false)
end

plot!(p1, t_full, pR_fit,
      lw=2.5, color=:orange, ls=:dash,
      label="fit médian  γ=$(round(γ_med, sigdigits=4))")

scatter!(p1, times, true_pR,
         ms=8, color=:tomato, markershape=:diamond,
         label="vraie aux timepoints")

p2 = histogram(gamma_errors,
               bins=15, color=:mediumseagreen, alpha=0.75,
               xlabel="Erreur relative γ (%)", ylabel="Compte",
               title="Distribution erreur γ  (n_mc=$(N_MC))",
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
          plot_title="Design : T=$(Int(T_MAX)) j | N=$(N_CELLS) | $(N_GROUPS) timepoints | $(N_REP) réplicats",
          plot_titlefontsize=20,
          guidefontsize=16, tickfontsize=14, legendfontsize=13)

savefig(pl, joinpath(@__DIR__, "gamma_single_design.png"))
display(pl)

println("\nFigure sauvegardée : gamma_single_design.png")