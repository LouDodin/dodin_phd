## ===== Packages =====
using DifferentialEquations
using OrdinaryDiffEqRosenbrock
using CSV
using DataFrames
using Statistics
using Plots
using Measures
using BlackBoxOptim
using LaTeXStrings
using Random
using SciMLBase
using Dates
using Colors
using DataInterpolations


## ===== Choices =====
const MODEL_FILE  = joinpath(@__DIR__, "models/SV_r.jl")
const replicates  = ["A", "B", "C"]
const cycles_fit  = 5
const cycles_sim  = 5
const n_runs      = 1
const n_interior  = [0,0,0,0,0]   # interior knots per cycle for c and p splines


## ===== Utils =====
include(MODEL_FILE)
using .ModelDef

output_dir = joinpath(@__DIR__, "output/r/$(n_interior[1])_$(n_interior[2])_$(n_interior[3])_$(n_interior[4])_$(n_interior[5])")
mkpath(output_dir)

const isoutofdomain = (u, p, t) -> any(x -> x < 0, u)

replicate_colors = Dict(
    "A" => RGB(0.6, 0.8, 1.0),
    "B" => RGB(31/255, 119/255, 180/255),
    "C" => RGB(0.0, 0.3, 0.7)
)
model_color = RGB(255/255, 127/255, 14/255)
cycle_color = RGB(31/255, 119/255, 180/255)


## ===== Input: S & V data =====

raw_data = Dict{String, Vector{NamedTuple}}()

for rep in replicates
    entries    = NamedTuple[]
    t_offset   = nothing
    t_end_prev = nothing

    for cyc_idx in 1:cycles_sim
        path_H = "modelisation/input/xp_input_20/hostData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cyc_idx).csv"
        path_V = "modelisation/input/xp_input_20/virusData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cyc_idx).csv"
        (!isfile(path_H) || !isfile(path_V)) && continue

        df_H = CSV.read(path_H, DataFrame)
        df_V = CSV.read(path_V, DataFrame)

        tH = df_H[:, 1] ./ 24.0;  H = Vector{Float64}(df_H[:, 2])
        tV = df_V[:, 1] ./ 24.0;  V = Vector{Float64}(df_V[:, 2])

        if t_offset === nothing; t_offset = tH[1]; end
        tH .-= t_offset;  tV .-= t_offset

        if t_end_prev !== nothing
            gap = tH[1] - t_end_prev
            tH .-= gap;  tV .-= gap
        end
        t_end_prev = tH[end]

        push!(entries, (index=cyc_idx, tH=tH, H=H, tV=tV, V=V))
    end
    raw_data[rep] = entries
end


## ===== Spline knot construction =====
#
# Knots are shared across all cycles and connected:
#   - t_start of cycle 1  (boundary)
#   - n_interior[1] equally-spaced interior knots inside cycle 1
#   - t_end of cycle 1 = t_start of cycle 2  (shared boundary)
#   - n_interior[2] interior knots inside cycle 2
#   - ...
#   - t_end of last cycle  (boundary)
#
# Total knots = 1 + sum(n_interior) + cycles_fit
# Each "inter-cycle" boundary node is shared (continuity).

function build_spline_knots(raw_data, replicates, cycles_fit, n_interior)
    knots = Float64[]

    for cyc in 1:cycles_fit
        t0 = minimum(raw_data[rep][cyc].tH[1] for rep in replicates)
        t1 = maximum(max(raw_data[rep][cyc].tH[end],
                         raw_data[rep][cyc].tV[end]) for rep in replicates)

        push!(knots, t0)   # left boundary of this cycle (= right boundary of previous)
        for j in 1:n_interior[cyc]
            push!(knots, t0 + j * (t1 - t0) / (n_interior[cyc] + 1))
        end
    end

    # Final right boundary (end of last cycle)
    cyc = cycles_fit
    t1_last = maximum(max(raw_data[rep][cyc].tH[end],
                          raw_data[rep][cyc].tV[end]) for rep in replicates)
    push!(knots, t1_last)

    return knots
end

const spline_knots = build_spline_knots(raw_data, replicates, cycles_fit, n_interior)
const N_KNOTS      = length(spline_knots)   # = 1 + sum(n_interior) + cycles_fit

println("Spline knots ($(N_KNOTS) total): ", round.(spline_knots, digits=2))


## ===== Spline builder =====
#
# c values are in log-space → exp to natural space
# p values are in logit-space → sigmoid to stay in (0,1)
#
# θ layout:
#   θ[1]           : log(ϕ)             (scalar)
#   θ[2:N_KNOTS+1] : log(c) at knots
#   θ[N_KNOTS+2 : 2*N_KNOTS+1] : logit(p) at knots

sigmoid(x) = 1.0 / (1.0 + exp(-x))

function build_spline_func(knots::Vector{Float64}, values::Vector{Float64})
    # values sont déjà dans l'espace de fitting (log pour c, logit pour p)
    # on interpole directement dans cet espace pour éviter les oscillations négatives
    itp = CubicSpline(values, knots)
    return t -> itp(clamp(t, knots[1], knots[end]))
end

function decode_theta(θ::Vector{Float64})
    ϕ = exp(θ[1])
    β = exp(θ[2])
    log_r = θ[3 : N_KNOTS+2]

    # Splines construites en espace transformé (log / logit)
    log_r_func   = build_spline_func(spline_knots, log_r)

    # Fonctions finales : transformation appliquée APRÈS interpolation
    r_func = t -> exp(log_r_func(t))

    # Valeurs aux nœuds pour le plot
    r_vals = exp.(log_r)

    return ϕ, β, r_func, r_vals
end


## ===== Search bounds =====
#
# From ModelDef.FITTED_PARAMS:
#   :ϕ  → log-space bounds for θ[1]
#   :c  → log-space bounds for θ[2..N_KNOTS+1]
#   :p  → logit-space bounds for θ[N_KNOTS+2..2*N_KNOTS+1]

fp_dict = Dict(fp.name => fp for fp in ModelDef.FITTED_PARAMS)

log_ϕ_lo  = log(fp_dict[:ϕ].lower);   log_ϕ_hi  = log(fp_dict[:ϕ].upper)
log_β_lo  = log(fp_dict[:β].lower);   log_β_hi  = log(fp_dict[:β].upper)
log_r_lo  = log(fp_dict[:r].lower);   log_r_hi  = log(fp_dict[:r].upper)

search_range = vcat(
    [(log_ϕ_lo, log_ϕ_hi)],
    [(log_β_lo, log_β_hi)],
    fill((log_r_lo, log_r_hi), N_KNOTS)
)

println("Total parameters: $(N_KNOTS+2)")


## ===== ODE solver =====

function solve_cycle(model!, ϕ, β, r_func, u0, t0, t1, t_save)
    params = (ϕ, β, r_func)
    prob   = ODEProblem(model!, u0, (t0, t1), params)
    sol    = solve(prob,
                   Rodas5(),
                   reltol = 1e-6,
                   abstol = 1e-6,
                   saveat = t_save,
                   isoutofdomain = isoutofdomain)
    return sol
end


## ===== Objective function =====

function objective(θ)
    ϕ, β, r_func, _ = decode_theta(θ)

    total_err = 0.0

    for rep in replicates
        for cycle in 1:cycles_fit
            data = raw_data[rep][cycle]
            H0   = data.H[1];  V0 = data.V[1]
            u0   = [H0, V0]

            t_data = sort(unique(vcat(data.tH, data.tV)))
            t0_data, t1_data = t_data[1], t_data[end]

            sol = solve_cycle(ModelDef.ODE_MODEL!, ϕ, β, r_func, u0, t0_data, t1_data, t_data)

            if sol.retcode != SciMLBase.ReturnCode.Success ||
               any(x -> x < 0, reduce(vcat, sol.u))
                return 1e15
            end

            H_pred = [sol(t)[1] for t in data.tH]
            V_pred = [sol(t)[2] for t in data.tV]

            total_err += sum((log.(H_pred) .- log.(data.H)).^2) / length(data.H)
            total_err += sum((log.(V_pred) .- log.(data.V)).^2) / length(data.V)
        end
    end

    return total_err
end


## ===== Optimisation =====

function run_DE(seed)
    Random.seed!(seed)
    res = bboptimize(
        objective;
        SearchRange          = search_range,
        NumDimensions        = N_KNOTS+2,
        Method               = :adaptive_de_rand_1_bin_radiuslimited,
        PopulationSize       = 250,
        MaxSteps             = 200_000,
        DifferentialWeight   = 0.9,
        crossoverProbability = 0.9,
        TraceMode            = :compact,
    )
    return (fitness = best_fitness(res), θ = best_candidate(res))
end

function main()
    results = Vector{NamedTuple}(undef, n_runs)
    Threads.@threads for i in 1:n_runs
        println("  Thread $(Threads.threadid()) — run $i (seed=$(1000+i))")
        results[i] = Base.invokelatest(run_DE, 1000 + i)
    end
    return results
end

results = Base.invokelatest(main)
best_idx  = argmin(res.fitness for res in results)
best_result = results[best_idx]
θbest     = best_result.θ

ϕ_best, β_best, r_func_best, r_knot_vals = decode_theta(θbest)

println("\n===== Best result ($(ModelDef.MODEL_NAME)) =====")
println("  Fitness = ", best_result.fitness)
println("  ϕ       = ", ϕ_best)
println("  β       = ", β_best)
println("  r at knots: ", round.(r_knot_vals, sigdigits=3))


## ===== Information Criteria =====

n_obs = sum(
    length(raw_data[rep][cyc].H) + length(raw_data[rep][cyc].V)
    for rep in replicates
    for cyc in 1:cycles_fit
)

function raw_rss(θ)
    ϕ, β, r_func, _ = decode_theta(θ)
    rss = 0.0

    for rep in replicates
        for cycle in 1:cycles_fit
            data = raw_data[rep][cycle]
            H0   = data.H[1];  V0 = data.V[1]
            u0   = [H0, V0]

            t_data = sort(unique(vcat(data.tH, data.tV)))
            sol    = solve_cycle(ModelDef.ODE_MODEL!, ϕ, β, r_func, u0, t_data[1], t_data[end], t_data)

            (sol.retcode != SciMLBase.ReturnCode.Success ||
             any(x -> x < 0, reduce(vcat, sol.u))) && return Inf

            H_pred = [sol(t)[1] for t in data.tH]
            V_pred = [sol(t)[2] for t in data.tV]

            rss += sum((log.(H_pred) .- log.(data.H)).^2)
            rss += sum((log.(V_pred) .- log.(data.V)).^2)
        end
    end
    return rss
end

rss_best = raw_rss(θbest)
k        = N_KNOTS+1
σ²       = rss_best / n_obs
logL     = -n_obs/2 * log(2π) - n_obs/2 * log(σ²) - n_obs/2
AIC      = 2k - 2logL
AICc     = AIC + 2k*(k+1) / (n_obs - k - 1)
BIC      = k*log(n_obs) - 2logL

println("\n===== Information Criteria ($(ModelDef.MODEL_NAME)) =====")
println("  n_obs = $n_obs,  k = $k")
println("  RSS   = $(round(rss_best, sigdigits=4))")
println("  logL  = $(round(logL,     sigdigits=6))")
println("  AIC   = $(round(AIC,      sigdigits=6))")
println("  AICc  = $(round(AICc,     sigdigits=6))")
println("  BIC   = $(round(BIC,      sigdigits=6))")


## ===== Simulation & Plot =====

ytick_vals_H   = [10.0^i for i in 2:9]
ytick_labels_H = [L"10^{%$i}" for i in 2:9]
ytick_vals_V   = [10.0^i for i in 5:10]
ytick_labels_V = [L"10^{%$i}" for i in 5:10]

title_str = "$(ModelDef.MODEL_NAME) - r(t) - ϕ=$(ϕ_best) - β=$(β_best)\n" *
            "Fitness=$(round(best_result.fitness, sigdigits=4))  " *
            "AIC=$(round(AIC, sigdigits=4))  AICc=$(round(AICc, sigdigits=4))  BIC=$(round(BIC, sigdigits=4))"

common_kw = (
    left_margin   = 15mm,
    right_margin  = 10mm,
    bottom_margin = 5mm,
    grid          = true,
    xlims         = (0,67),
    xtickfontsize = 20,
    ytickfontsize = 20,
    legendfontsize= 20,
    guidefontsize = 20,
    xlabel        = "Time [days]",
    legend        = :bottomright,
)

pl_fit = plot(
    layout            = (3, 1),
    size              = (2500, 2000),
    top_margin        = 10mm,
    plot_title        = title_str,
    plot_titlefontsize= 25,
)

for rep in replicates
    for cyc in 1:cycles_sim
        data = raw_data[rep][cyc]
        scatter!(pl_fit[1], data.tH, data.H;
            color = replicate_colors[rep],
            label = cyc == 1 ? "Rep $rep" : "",
            ylabel = "Host [cell/mL]",
            yscale = :log10, ylims = (1e2, 1e9),
            yticks = (ytick_vals_H, ytick_labels_H),
            markersize = 10, alpha = 0.7,
            common_kw..., legend=(0.15, 0.9)
        )
    end
end

for rep in replicates
    for cyc in 1:cycles_sim
        data = raw_data[rep][cyc]
        scatter!(pl_fit[2], data.tV, data.V;
            color = replicate_colors[rep],
            label = cyc == 1 ? "Rep $rep" : "",
            ylabel = "Virus [part/mL]",
            yscale = :log10, ylims = (1e5, 1e10),
            yticks = (ytick_vals_V, ytick_labels_V),
            markersize = 10, alpha = 0.7,
            common_kw..., legend=(0.15, 0.3)
        )
    end
end

for cycle in 1:cycles_sim
    H0 = mean(raw_data[rep][cycle].H[1] for rep in replicates)
    V0 = mean(raw_data[rep][cycle].V[1] for rep in replicates)
    u0 = [H0, V0]

    t0 = minimum(raw_data[rep][cycle].tH[1]   for rep in replicates)
    t1 = maximum(max(raw_data[rep][cycle].tH[end],
                     raw_data[rep][cycle].tV[end]) for rep in replicates)

    sol = solve(
        ODEProblem(ModelDef.ODE_MODEL!, u0, (t0, t1), (ϕ_best, β_best, r_func_best)),
        Rodas5(), reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain,
    )

    S_traj = [u[1] for u in sol.u]
    V_traj = [u[2] for u in sol.u]
    r_traj = r_func_best.(sol.t)

    lbl1 = cycle == 1
    plot!(pl_fit[1], sol.t, S_traj; label=lbl1 ? "H model" : "", color=model_color, lw=8, alpha=0.8)
    plot!(pl_fit[2], sol.t, V_traj; label=lbl1 ? "V model" : "", color=model_color,  lw=8, alpha=0.8)
    plot!(pl_fit[3], sol.t, r_traj; label=lbl1 ? "r model" : "", color=model_color,  lw=8, alpha=0.8)
end

scatter!(pl_fit[3], spline_knots, r_knot_vals;
    color = :red, markersize = 10, markershape = :circle, label = "r knots",
            common_kw..., ylabel="r\n[/day]"
)

# ── Dilution lines on all 3 panels ───────────────────────────────────────
cycle_tbounds = [(
    minimum(raw_data[rep][cyc].tH[1] for rep in replicates),
    maximum(max(raw_data[rep][cyc].tH[end], raw_data[rep][cyc].tV[end]) for rep in replicates)
) for cyc in 1:cycles_sim]

dilution_times = [cycle_tbounds[c][2] for c in 1:cycles_sim-1]

for panel in 1:3
    for td in dilution_times
        vline!(pl_fit[panel], [td]; color=cycle_color, lw=4, ls=:dot, label=nothing)
    end
end

display(pl_fit)

fig_path = "$(output_dir)/$(ModelDef.MODEL_NAME)_spline_plot.png"
savefig(pl_fit, fig_path)
println("\nFigure saved to $fig_path")


## ===== Log =====
log_path = "$(output_dir)/$(ModelDef.MODEL_NAME)_spline_log.txt"
open(log_path, "w") do io
    println(io, "===== Run log =====")
    println(io, "Date        : ", Dates.now())
    println(io, "Model       : ", ModelDef.MODEL_NAME)
    println(io, "Model file  : ", MODEL_FILE)
    println(io, "n_runs      : ", n_runs)
    println(io, "cycles_fit  : ", cycles_fit)
    println(io, "cycles_sim  : ", cycles_sim)
    println(io, "replicates  : ", replicates)
    println(io, "n_interior  : ", n_interior)
    println(io, "N_KNOTS     : ", N_KNOTS)
    println(io)
    println(io, "===== Fixed parameters =====")
    for (k, v) in pairs(ModelDef.FIXED_PARAMS)
        println(io, "  $k = $v")
    end
    println(io)
    println(io, "===== Spline knots =====")
    for (i, tk) in enumerate(spline_knots)
        println(io, "  knot[$i] = $(round(tk, digits=4)) days")
    end
    println(io)
    println(io, "===== All runs =====")
    for (i, res) in enumerate(results)
        println(io, "  Run $i (seed=$(1000+i)) : loss=$(res.fitness)")
    end
    println(io)
    println(io, "===== Best result =====")
    println(io, "  Loss = ", best_result.fitness)
    println(io, "  ϕ         : ", ϕ_best)
    println(io, "  β         : ", β_best)
    println(io, "  r at knots : ", r_knot_vals)
    println(io)
    println(io, "===== Information Criteria =====")
    println(io, "  n_obs = $n_obs")
    println(io, "  k     = $k  (fitted params)")
    println(io, "  RSS   = $rss_best")
    println(io, "  logL  = $logL")
    println(io, "  AIC   = $AIC")
    println(io, "  AICc  = $AICc")
    println(io, "  BIC   = $BIC")
end
println("Log saved to $log_path")