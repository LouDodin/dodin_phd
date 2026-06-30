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

# ── Simulation ODE — grille incluant exactement les timepoints ────────────────
# query_times : vecteur des temps expérimentaux à retourner exactement
# retourne (t_grid, pR_grid, pR_at_query) 
function simulate(t_end, gamma, query_times; dt=0.5)
    # Grille temporelle = pas régulier + timepoints expérimentaux
    t_grid = sort(unique(vcat(0.0, 0.0:dt:t_end, query_times)))

    S = (1 - prop_R0) * H0
    R = prop_R0 * H0

    pR_grid  = zeros(length(t_grid))
    pR_grid[1] = R / (S + R)

    for i in 2:length(t_grid)
        dt_i = t_grid[i] - t_grid[i-1]   # pas variable selon la grille
        H     = S + R
        S_new = max(0.0, S + dt_i * (TRUE_r * S * (1 - H/TRUE_K) - TRUE_alpha * S + gamma * R))
        R_new = max(0.0, R + dt_i * (TRUE_r * R * (1 - H/TRUE_K) + TRUE_alpha * S - gamma * R))
        S, R  = S_new, R_new
        pR_grid[i] = R / (S + R)
    end

    # Extraire directement les valeurs aux timepoints expérimentaux
    query_idx  = [findfirst(==(t), t_grid) for t in query_times]
    pR_at_query = pR_grid[query_idx]

    return t_grid, pR_grid, pR_at_query
end

# ── Échantillonnage binomial (mesure par isolation clonale) ───────────────────
sample_binomial(n, p) = sum(rand() < p for _ in 1:n) / n

# ── Fitting — 1 paramètre γ, BlackBoxOptim DE ────────────────────────────────
function fit_gamma(times, obs)
    function objective(x)
        gamma = x[1]
        _, _, pR_pred = simulate(maximum(times), gamma, times)
        return sum((obs .- pR_pred).^2)
    end

    res = bboptimize(
        objective;
        SearchRange          = [(1e-4, 0.5)],   # bornes sur γ
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

# ── Monte Carlo ───────────────────────────────────────────────────────────────
function run_monte_carlo(times, true_pR; n_cells, n_rep, n_mc)
    gamma_ests   = Float64[]
    gamma_errors = Float64[]

    for mc in 1:n_mc
        # Moyenne de n_rep réplicats binomiaux à chaque timepoint
        obs = [mean(sample_binomial(n_cells, true_pR[i]) for _ in 1:n_rep)
               for i in eachindex(times)]

        γ_est = fit_gamma(times, obs)
        push!(gamma_ests,   γ_est)
        push!(gamma_errors, abs(γ_est - TRUE_gamma) / TRUE_gamma * 100)

        @printf("  MC %2d/%d — γ_est=%.5f  erreur=%.1f%%\n",
                mc, n_mc, γ_est, gamma_errors[end])
    end

    return gamma_ests, gamma_errors
end

# ── Plots ─────────────────────────────────────────────────────────────────────
function make_plots(; n_groups=8, t_max=150.0, n_cells=100, n_rep=3, n_mc=50, seed=42)
    Random.seed!(seed)

    times = [t_max * k / n_groups for k in 1:n_groups]

    # Valeurs vraies exactement aux timepoints (pas d'interpolation)
    _, _, true_pR = simulate(t_max, TRUE_gamma, times)

    println("=== Estimation de γ ===")
    @printf("n_groups=%d | t_max=%.0f j | n_cells=%d | n_rep=%d | n_mc=%d\n\n",
            n_groups, t_max, n_cells, n_rep, n_mc)

    gamma_ests, gamma_errors = run_monte_carlo(times, true_pR; n_cells, n_rep, n_mc)

    γ_med = median(gamma_ests)
    err_med = median(gamma_errors)

    @printf("\nRésultats (médiane sur %d MC) :\n", n_mc)
    @printf("  γ vrai   = %.6f\n", TRUE_gamma)
    @printf("  γ estimé = %.6f\n", γ_med)
    @printf("  Erreur   = %.1f%%\n", err_med)

    # Courbes pour le plot
    t_full, pR_full, _ = simulate(t_max, TRUE_gamma, times)
    _, pR_fit, _        = simulate(t_max, γ_med,     times)

    # ── P1 : courbe vraie + fit + timepoints ──────────────────────────────
    p1 = plot(t_full, pR_full,
              lw=3, color=:steelblue,
              label="vraie (γ=$(TRUE_gamma))",
              xlabel="Temps [jours]", ylabel="prop_R [-]",
              title="prop_R(t)", ylims=(0, 1.05), legend=:bottomleft)
    plot!(p1, t_full, pR_fit,
          lw=2, color=:orange, ls=:dash,
          label="fit médian (γ=$(round(γ_med, sigdigits=4)))")
    scatter!(p1, times, true_pR,
             ms=7, color=:tomato, label="timepoints")
    eq = TRUE_alpha / (TRUE_alpha + TRUE_gamma)
    hline!(p1, [eq], lw=1, ls=:dot, color=:gray,
           label="équilibre=$(round(eq, sigdigits=3))")

    # ── P2 : histogramme erreur γ ──────────────────────────────────────────
    p2 = histogram(gamma_errors,
                   bins=20, color=:mediumseagreen, alpha=0.75,
                   xlabel="Erreur relative γ (%)", ylabel="Compte",
                   title="Distribution erreur γ  (n_mc=$n_mc)", label=false)
    vline!(p2, [err_med], lw=2, color=:white,
           label="médiane=$(round(err_med, digits=1))%")

    # ── P3 : erreur γ vs N cellules ────────────────────────────────────────
    n_vals    = [10, 20, 50, 100, 200, 500]
    med_errs  = Float64[]
    p25_errs  = Float64[]
    p75_errs  = Float64[]

    println("\nScan N cellules...")
    for n in n_vals
        _, errs = run_monte_carlo(times, true_pR; n_cells=n, n_rep, n_mc)
        push!(med_errs, median(errs))
        push!(p25_errs, quantile(errs, 0.25))
        push!(p75_errs, quantile(errs, 0.75))
        @printf("  N=%3d → erreur γ médiane = %.1f%%\n", n, med_errs[end])
    end

    p3 = plot(n_vals, med_errs,
              lw=3, marker=:circle, ms=6, color=:mediumseagreen,
              xlabel="N cellules par puit", ylabel="Erreur γ (%)",
              title="Erreur γ vs N cellules", label="médiane", legend=:topright)
    plot!(p3, n_vals, p25_errs, fillrange=p75_errs,
          alpha=0.2, color=:mediumseagreen, label="IQR")
    hline!(p3, [20.0], lw=1, ls=:dash, color=:red,    label="seuil 20%")
    hline!(p3, [10.0], lw=1, ls=:dash, color=:orange, label="seuil 10%")

    # ── P4 : erreur γ vs T_max ────────────────────────────────────────────
    t_vals   = [30, 50, 75, 100, 150, 200]
    err_tmax = Float64[]

    println("\nScan T_max...")
    for t in t_vals
        ts = [t * k / n_groups for k in 1:n_groups]
        _, _, tv = simulate(t, TRUE_gamma, ts)
        _, errs  = run_monte_carlo(ts, tv; n_cells, n_rep, n_mc)
        push!(err_tmax, median(errs))
        @printf("  T=%3d j → erreur γ médiane = %.1f%%\n", t, err_tmax[end])
    end

    p4 = plot(t_vals, err_tmax,
              lw=3, marker=:circle, ms=6, color=:cornflowerblue,
              xlabel="Temps total [jours]", ylabel="Erreur γ (%)",
              title="Erreur γ vs Temps total", label="médiane", legend=:topright)
    hline!(p4, [20.0], lw=1, ls=:dash, color=:red,    label="seuil 20%")
    hline!(p4, [10.0], lw=1, ls=:dash, color=:orange, label="seuil 10%")

    pl = plot(p1, p2, p3, p4,
              layout=(2, 2),
              size=(2400, 1600),
              left_margin=15mm, right_margin=8mm,
              top_margin=8mm,   bottom_margin=12mm,
              plot_title="Estimation de γ — Design expérimental",
              plot_titlefontsize=22, titlefontsize=16,
              guidefontsize=14, tickfontsize=12, legendfontsize=11)

    savefig(pl, joinpath(@__DIR__, "gamma_design.png"))
    display(pl)
    println("\nFigure sauvegardée : gamma_design.png")
    return pl
end

# ── Point d'entrée ────────────────────────────────────────────────────────────
make_plots(
    n_groups = 8,
    t_max    = 150.0,
    n_cells  = 100,
    n_rep    = 3,
    n_mc     = 50,
)