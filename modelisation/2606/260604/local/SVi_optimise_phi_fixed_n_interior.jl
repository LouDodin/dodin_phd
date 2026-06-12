## ===== Import packages =====
using Dates
using CSV
using DataFrames
using DifferentialEquations
using SciMLBase
using Statistics
using Plots
using Measures
using BlackBoxOptim
using Base.Threads
using DataInterpolations
using Printf
using LaTeXStrings
using OrdinaryDiffEqRosenbrock


## ===== CONFIGURATION =====

const DATA_DIR   = joinpath(@__DIR__, "../../../input/xp_input_20")
const OUTPUT_DIR = joinpath(@__DIR__, "output")
const REPLICATES = ["A", "B", "C"]
const N_RUNS     = 1
const N_INTERIOR = [0, 0, 0, 0, 0]
const N_CYCLES   = length(N_INTERIOR)

# Paramètres biologiques (fixés)
const r = 0.574619342477644
const K = 6.675449070379925e7
const β = 144.0
const δ = 0.02

# Couleurs
const COLOR_REPS  = [RGB(0.6, 0.8, 1.0), RGB(31/255, 119/255, 180/255), RGB(0.0, 0.3, 0.7)]
const COLOR_MODEL = RGB(255/255, 127/255, 14/255)

yticks_log(lo, hi) = (
    [10.0^i for i in lo:hi],
    [L"10^{%$i}" for i in lo:hi]
)


## ===== STRUCTURE =====

struct Cycle
    rep   :: String
    index :: Int
    tH    :: Vector{Float64}
    H     :: Vector{Float64}
    tV    :: Vector{Float64}
    V     :: Vector{Float64}
    u0    :: Vector{Float64}
end


## ===== IMPORT =====

function load_replicate(rep::String)::Vector{Cycle}
    cycles     = Cycle[]
    t_offset   = nothing   # ancrage global à t=0 (cycle 1)
    t_end_prev = nothing   # fin du cycle précédent (pour coller les cycles)

    for cyc in 1:N_CYCLES
        df_H = CSV.read(joinpath(DATA_DIR,
            "hostData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cyc).csv"),
            DataFrame)
        df_V = CSV.read(joinpath(DATA_DIR,
            "virusData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cyc).csv"),
            DataFrame)

        tH = df_H[:, 1] ./ 24.0
        H  = df_H[:, 2]
        tV = df_V[:, 1] ./ 24.0
        V  = df_V[:, 2]

        # 1) Ancrer le cycle 1 à t = 0
        if t_offset === nothing
            t_offset = tH[1]
        end
        tH .-= t_offset
        tV .-= t_offset

        # 2) Coller les cycles bout à bout (supprimer le gap inter-cycles)
        if t_end_prev !== nothing
            gap = tH[1] - t_end_prev
            tH .-= gap
            tV .-= gap
        end
        t_end_prev = tH[end]

        push!(cycles, Cycle(rep, cyc, tH, H, tV, V, [H[1], V[1]]))
    end
    return cycles
end

all_cycles     = vcat([load_replicate(rep) for rep in REPLICATES]...)
cycles_per_rep = Dict(rep => filter(c -> c.rep == rep, all_cycles) for rep in REPLICATES)

println("Cycles chargés : $(length(all_cycles)) au total")
for rep in REPLICATES
    cs = cycles_per_rep[rep]
    println("  Rep $rep : t ∈ [$(round(cs[1].tH[1], digits=2)), $(round(cs[end].tH[end], digits=2))] jours")
end


## ===== CONSTRUCTION DES NŒUDS =====

function build_global_knots()::Vector{Float64}
    knots = Float64[]
    for cyc in 1:N_CYCLES
        cyc_ref = cycles_per_rep["A"][cyc]
        t0 = cyc_ref.tH[1]
        t1 = cyc_ref.tH[end]

        push!(knots, t0)
        n_int = N_INTERIOR[cyc]
        if n_int > 0
            inner = collect(range(t0, t1, length=n_int+2)[2:end-1])
            append!(knots, inner)
        end
        push!(knots, t1)
    end
    return sort(unique(knots))
end

global_knots = build_global_knots()
n_knots      = length(global_knots)

println("\nNœuds globaux ($n_knots au total) :")
println("  ", round.(global_knots, digits=2))
println("  N_INTERIOR par cycle : $N_INTERIOR")
println("  Dimension de θ : $n_knots\n")


## ===== MODÈLE ODE =====

function ode_model!(dY, Y, phi_func, t)
    S, Vi = Y[1], Y[2]
    ϕ = phi_func(t)
    dY[1] = r * S * (1 - S/K) - ϕ * S * Vi
    dY[2] = β * ϕ * S * Vi - δ * Vi
end

function integrate_cycle(cyc::Cycle, phi_func)
    t0 = cyc.tH[1]
    t1 = max(cyc.tH[end], cyc.tV[end])
    sv = sort(unique(vcat(cyc.tH, cyc.tV)))

    prob = ODEProblem(ode_model!, cyc.u0, (t0, t1), phi_func)
    sol  = solve(prob, Rodas5(),
                 reltol=1e-6, abstol=1e-6,
                 saveat=sv,
                 isoutofdomain=(u, p, t) -> any(x -> x < 0, u))

    return SciMLBase.successful_retcode(sol) ? sol : nothing
end


## ===== FONCTION OBJECTIF =====

function objective(θ::Vector{Float64}, t_knots::Vector{Float64})
    spline   = CubicSpline(θ, t_knots)
    phi_func = t -> exp(spline(clamp(t, t_knots[1], t_knots[end])))
    total_err = 0.0

    for cyc in all_cycles
        sol = integrate_cycle(cyc, phi_func)
        sol === nothing && return 1e12

        S_pred = [max(sol(t)[1], 1e-12) for t in cyc.tH]
        V_pred = [max(sol(t)[2], 1e-12) for t in cyc.tV]

        total_err += sum((log.(S_pred) .- log.(cyc.H)).^2) / length(cyc.H)
        total_err += sum((log.(V_pred) .- log.(cyc.V)).^2) / length(cyc.V)
    end
    return total_err
end


## ===== OPTIMISATION =====

function run_DE(seed::Int, t_knots::Vector{Float64})
    n  = length(t_knots)
    lo = fill(log(1e-15), n)
    hi = fill(log(1e-6),  n)

    res = bboptimize(
        θ -> objective(θ, t_knots);
        SearchRange          = collect(zip(lo, hi)),
        NumDimensions        = n,
        Method               = :xnes,
        PopulationSize       = 1000,
        MaxSteps             = 10000,
        DifferentialWeight   = 0.5,
        CrossoverProbability = 0.9,
        TraceMode            = :silent,
        RandomSeed           = seed,
    )
    return (fitness=best_fitness(res), θ=best_candidate(res), seed=seed)
end


## ===== GRAPHIQUE =====

function make_global_plot(t_knots, θbest, best_fit, best_seed)
    spline   = CubicSpline(θbest, t_knots)
    phi_func = t -> exp(spline(clamp(t, t_knots[1], t_knots[end])))

    t_fine    = range(t_knots[1], t_knots[end], length=2000)
    phi_curve = phi_func.(t_fine)

    yt_host  = yticks_log(3, 9)
    yt_virus = yticks_log(6, 10)
    exp_lo   = floor(Int, log10(minimum(phi_curve)))
    exp_hi   = ceil(Int,  log10(maximum(phi_curve)))
    phi_lo   = 10.0^exp_lo
    phi_hi   = 10.0^exp_hi

    t_lo = t_knots[1]  - 1
    t_hi = t_knots[end] + 1

    pl = plot(
        layout=(3, 1), size=(1600, 1200),
        left_margin=18mm, right_margin=10mm,
        top_margin=15mm,  bottom_margin=10mm,
        grid=true,
        xlims=(t_lo, t_hi),
        ytickfontsize=18, legendfontsize=14,
        guidefontsize=18, xtickfontsize=16,
        titlefontsize=13, xlabel="Time (days)",
        legend=:bottomright,
    )

    title_str = "Global fit — fitness=$(round(best_fit, sigdigits=5)) | seed=$best_seed | n_interior=$(N_INTERIOR)"
    plot!(pl, title=title_str, subplot=1)

    # Bandes de fond par cycle
    for cyc_idx in 1:N_CYCLES
        cyc_ref = cycles_per_rep["A"][cyc_idx]
        t0c = cyc_ref.tH[1]
        t1c = cyc_ref.tH[end]
        fc  = cyc_idx % 2 == 0 ? RGB(0.93, 0.93, 0.93) : :white
        for sp in 1:3
            vspan!(pl[sp], [t0c, t1c]; color=fc, alpha=0.35, label=false)
        end
        annotate!(pl[3], (t0c + t1c) / 2, phi_hi * 0.8,
                  text("C$cyc_idx", 11, :center, :gray))
    end

    # Données des 3 réplicats
    for (i, rep) in enumerate(REPLICATES)
        tH = vcat([c.tH for c in cycles_per_rep[rep]]...)
        H  = vcat([c.H  for c in cycles_per_rep[rep]]...)
        tV = vcat([c.tV for c in cycles_per_rep[rep]]...)
        V  = vcat([c.V  for c in cycles_per_rep[rep]]...)

        scatter!(pl[1], tH, H;
                 color=COLOR_REPS[i], alpha=0.7, label="Rep $rep",
                 ylabel="Host abundance (cell/mL)",
                 ylims=(1e3, 1e9), yticks=yt_host, yscale=:log10,
                 markershape=:circle, markersize=6, markerstrokewidth=1)
        scatter!(pl[2], tV, V;
                 color=COLOR_REPS[i], alpha=0.7, label="Rep $rep",
                 ylabel="Virus abundance (virion/mL)",
                 ylims=(1e6, 1e10), yticks=yt_virus, yscale=:log10,
                 markershape=:circle, markersize=6, markerstrokewidth=1)
    end

    # Trajectoires modèle (u0 = moyenne des 3 réplicats par cycle)
    for cyc_idx in 1:N_CYCLES
        cycs_this = filter(c -> c.index == cyc_idx, all_cycles)
        u0_mean   = vec(mean(hcat([c.u0 for c in cycs_this]...), dims=2))
        t0c = cycs_this[1].tH[1]
        t1c = maximum(max(c.tH[end], c.tV[end]) for c in cycs_this)

        prob = ODEProblem(ode_model!, u0_mean, (t0c, t1c), phi_func)
        sol  = solve(prob, Rodas5(), reltol=1e-6, abstol=1e-6,
                     isoutofdomain=(u, p, t) -> any(x -> x < 0, u))
        SciMLBase.successful_retcode(sol) || continue

        lbl = cyc_idx == 1 ? "Model" : false
        plot!(pl[1], sol.t, sol[1,:]; lw=2.5, color=COLOR_MODEL, label=lbl)
        plot!(pl[2], sol.t, sol[2,:]; lw=2.5, color=COLOR_MODEL, label=lbl)
    end

    # φ(t)
    scatter!(pl[3], t_knots, exp.(θbest);
             label="Knots", markershape=:circle, markersize=8,
             markerstrokewidth=0, color=COLOR_MODEL)
    plot!(pl[3], collect(t_fine), phi_curve;
          label="φ(t)", lw=3, color=COLOR_MODEL,
          ylabel="φ (mL/(cell·day))",
          yscale=:log10,
          ylims=(phi_lo, phi_hi),
          yticks=(10.0 .^ (exp_lo:exp_hi),
                  [L"10^{%$i}" for i in exp_lo:exp_hi]))

    return pl
end


## ===== EXPORT POLYNÔME =====

function export_polynomial(path, spline, t_knots, θbest, best_fit, best_seed, all_results)
    T    = spline.t
    A_sp = spline.u
    z    = spline.z

    open(path, "w") do io
        @printf(io, "========== METADATA ==========\n")
        @printf(io, "export_date      = %s\n",   string(now()))
        @printf(io, "n_cycles         = %d\n",   N_CYCLES)
        @printf(io, "replicates       = %s (joint fit)\n", join(REPLICATES, ", "))
        @printf(io, "n_interior/cycle = %s\n",   string(N_INTERIOR))
        @printf(io, "total_knots      = %d\n",   length(t_knots))
        @printf(io, "n_runs           = %d\n",   N_RUNS)
        @printf(io, "best_fitness     = %.8e\n", best_fit)
        @printf(io, "best_seed        = %d\n\n", best_seed)

        @printf(io, "========== ALL RUNS ==========\n")
        @printf(io, "%-6s  %-14s  %s\n", "seed", "fitness", "θ (log-space)")
        for res in sort(collect(all_results), by=r -> r.fitness)
            θ_str = join([@sprintf("%.4e", v) for v in res.θ], "  ")
            @printf(io, "%-6d  %-14.8e  [%s]\n", res.seed, res.fitness, θ_str)
        end
        @printf(io, "\n")

        @printf(io, "========== OPTIMAL KNOTS ==========\n")
        @printf(io, "%-4s  %-14s  %-14s  %s\n", "i", "t_i (days)", "log(φ(t_i))", "φ(t_i) (mL/cell/day)")
        for (i, (t, lv)) in enumerate(zip(t_knots, θbest))
            @printf(io, "%-4d  %-14.4f  %-14.6e  %-14.6e\n", i, t, lv, exp(lv))
        end
        @printf(io, "\n")

        @printf(io, "========== PIECEWISE POLYNOMIAL OF log(φ(t)) ==========\n")
        @printf(io, "φ(t) = exp(P(t)),  P(t) spline cubique C² sur %d intervalles\n\n", length(T)-1)

        for i in 1:(length(T)-1)
            ti, ti1 = T[i], T[i+1]
            h = ti1 - ti
            a = A_sp[i]
            b = (A_sp[i+1] - A_sp[i])/h - h*(2*z[i] + z[i+1])/6
            c = z[i]/2
            d = (z[i+1] - z[i])/(6*h)

            @printf(io, "Intervalle [%.4f, %.4f] jours :\n", ti, ti1)
            @printf(io, "  log(φ(t)) = %.6e\n",                 a)
            @printf(io, "           + %.6e * (t - %.4f)\n",     b, ti)
            @printf(io, "           + %.6e * (t - %.4f)^2\n",   c, ti)
            @printf(io, "           + %.6e * (t - %.4f)^3\n\n", d, ti)
        end
    end
end


## ===== MAIN =====

mkpath(OUTPUT_DIR)
println("Threads disponibles : $(Threads.nthreads())")
println("Dimension de θ : $n_knots paramètres\n")

# Optimisation parallèle
results = Vector{NamedTuple{(:fitness, :θ, :seed), Tuple{Float64, Vector{Float64}, Int}}}(undef, N_RUNS)
Threads.@threads for i in 1:N_RUNS
    println("  Thread $(threadid()) — run $i (seed=$i)")
    results[i] = Base.invokelatest(run_DE, i, global_knots)
end

# Meilleur résultat
best_idx = argmin(r.fitness for r in results)
best     = results[best_idx]
println("\nBest fitness = $(round(best.fitness, sigdigits=6))  (seed=$(best.seed))")
println("Tous les runs :")
for res in sort(collect(results), by=r -> r.fitness)
    println("  seed=$(res.seed)  fitness=$(round(res.fitness, sigdigits=6))")
end

# Graphique
pl        = make_global_plot(global_knots, best.θ, best.fitness, best.seed)
plot_path = joinpath(OUTPUT_DIR, "global_nint$(join(N_INTERIOR, '_'))_plot.png")
savefig(pl, plot_path)
println("\nGraphique sauvegardé → $plot_path")

# Export polynôme1
spline    = CubicSpline(best.θ, global_knots)
poly_path = joinpath(OUTPUT_DIR, "global_nint$(join(N_INTERIOR, '_'))_polynomial.txt")
export_polynomial(poly_path, spline, global_knots, best.θ,
                  best.fitness, best.seed, results)
println("Polynôme sauvegardé → $poly_path")