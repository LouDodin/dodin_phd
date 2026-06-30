using Plots
using Measures
using LaTeXStrings
using Statistics
using Random
using Printf
using BlackBoxOptim
using Distributions
using Base.Threads

# ─────────────────────────────────────────────────────────────────────────────
# PARAMÈTRES FIXES
# ─────────────────────────────────────────────────────────────────────────────
const TRUE_r     = 0.574619342477644
const TRUE_K     = 6.675449070379925e7
const TRUE_alpha = 2.185e-5
const H0         = 1e5
const prop_R0    = 1.0 - 1e-15

const T_TOTAL   = 20
const N_TIME    = 1
const N_REP     = 5
const N_CELLS   = 96
const N_MC      = 300

Random.seed!(42)

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
# ÉCHANTILLONNAGE BINOMIAL
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
# BOUCLE SUR γ VRAI
# ─────────────────────────────────────────────────────────────────────────────
γ_true_list = vcat(
    [9.47815e-3],
    10 .^ range(log10(1e-5), log10(1e-1), length=50)
)

median_errors = zeros(length(γ_true_list))
gamma_median_predicted = zeros(length(γ_true_list))

times = [T_TOTAL/N_TIME*i for i in 1:N_TIME]

for (k, TRUE_gamma) in enumerate(γ_true_list)

    println("TRUE_gamma = ", TRUE_gamma)

    _, _, true_pR = simulate(T_TOTAL, TRUE_gamma, times)

    γ_est   = zeros(N_MC)
    γ_error = zeros(N_MC)

    rngs = [MersenneTwister(1234 + i) for i in 1:N_MC]

    @threads for i in 1:N_MC
        rng = rngs[i]

        obs = [
            mean(sample_binomial(N_CELLS, true_pR[j], rng) for _ in 1:N_REP)
            for j in eachindex(times)
        ]

        γ_est[i] = fit_gamma(times, obs)
        γ_error[i] = abs(γ_est[i] - TRUE_gamma) / TRUE_gamma * 100
    end

    median_errors[k] = median(γ_error)
    gamma_median_predicted[k] = median(γ_est)

    println("Erreur médiane = ", median_errors[k], "%")
    println("Gamma prédit médian = ", gamma_median_predicted[k])
end

# ─────────────────────────────────────────────────────────────────────────────
# SAUVEGARDE CSV
# ─────────────────────────────────────────────────────────────────────────────
csv_path = joinpath(@__DIR__, "median_error_vs_gamma.csv")

open(csv_path, "w") do io
    println(io, "gamma_true,median_error_percent,gamma_median_predicted")
    for i in eachindex(γ_true_list)
        println(io, "$(γ_true_list[i]),$(median_errors[i]),$(gamma_median_predicted[i])")
    end
end

println("\nCSV sauvegardé : ", csv_path)