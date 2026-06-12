## =============================================================================
## main.jl — Generic fitting & simulation driver
##
## Usage:
##   julia main.jl                          # uses MODEL_FILE defined below
##   julia main.jl models/model_SRVi.jl    # or pass model file as CLI arg
##
## Each model file must expose a `ModelDef` module with the interface described
## in models/model_SRVi.jl.
##
## Extended interface (optional fields — backward-compatible) :
##   ModelDef.LOG_INDICES  Vector{Int}  indices of θ_opt where exp() is applied
##                         during decode_params.  If absent, behaviour falls back
##                         to the old LOG_FIT flag (exp applied to all or none).
##   ModelDef.PHI_equiv(Y_t, p, t)   3-argument form required when φ depends on t
##                         (e.g. SRVi_phi).  The 2-argument form is still accepted
##                         for backward compatibility.
## =============================================================================

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


## ===== Select model =====
const MODEL_FILE = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "models/model_SRVi_phi.jl")

println("Loading model from: $MODEL_FILE")
include(MODEL_FILE)
using .ModelDef


## ===== Choices =====
const ϕ_combi     = "3_2_2_3_2"
const replicates  = ["A", "B", "C"]
const cycles_fit  = 5
const cycles_sim  = 5
const n_runs      = 1

const λ = 1


## ===== Inputs =====

# --- S & V data ---
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

# Check plot
replicate_colors = Dict(
    "A" => RGB(0.6, 0.8, 1.0),
    "B" => RGB(31/255, 119/255, 180/255),
    "C" => RGB(0.0, 0.3, 0.7)
)

pl_data = plot(layout=(1,2), size=(900,350), margins=5mm, legend=:bottomright)
for rep in replicates
    for cyc in 1:cycles_sim
        data = raw_data[rep][cyc]
        scatter!(pl_data[1], data.tH, data.H;
            color=replicate_colors[rep], label=cyc==1 ? "Rep $rep" : "",
            xlabel="Time (days)", ylabel="Host abundance (cell/mL)",
            yscale=:log10, ylims=(1e2,1e8), title="H")
        scatter!(pl_data[2], data.tV, data.V;
            color=replicate_colors[rep], label=cyc==1 ? "Rep $rep" : "",
            xlabel="Time (days)", ylabel="Virus abundance (part/mL)",
            yscale=:log10, ylims=(1e3,1e10), title="Vi")
    end
end
display(pl_data)


# --- φ'(t) polynomial ---
const poly_file = joinpath(@__DIR__, "../../260604/genotoul/output/nint_$(ϕ_combi)/polynomial.txt")

function parse_polynome(filepath::String)
    isfile(filepath) || error("Polynomial file not found: $filepath")
    lines = readlines(filepath)

    metadata = Dict{String,Any}()
    for line in lines
        m = match(r"^(n_knots|n_intervalles|best_fitness|best_seed|n_runs)\s*=\s*(.+)$", strip(line))
        m === nothing && continue
        key = m.captures[1]
        val = strip(m.captures[2])
        metadata[key] = key in ("n_knots","n_intervalles","best_seed","n_runs") ?
                         parse(Int, val) : parse(Float64, val)
    end

    function extract_float(s::AbstractString)::Float64
        m = match(r"([+-]?\s*[0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)", strip(s))
        m === nothing && error("Cannot extract float from: \"$s\"")
        parse(Float64, replace(m.captures[1], " " => ""))
    end

    intervals = Vector{NTuple{6,Float64}}()
    i = 1
    while i <= length(lines)
        m_iv = match(r"Interval\s+\[([0-9eE+\-.]+),\s*([0-9eE+\-.]+)\]\s+days:", strip(lines[i]))
        if m_iv !== nothing
            if i + 4 > length(lines)
                println("Incomplete interval block at line $i — skipping")
                i += 1; continue
            end
            t0 = parse(Float64, m_iv.captures[1])
            t1 = parse(Float64, m_iv.captures[2])
            a  = extract_float(split(lines[i+1], "=")[end])
            b  = extract_float(split(lines[i+2], "*")[1])
            c  = extract_float(split(lines[i+3], "*")[1])
            d  = extract_float(split(lines[i+4], "*")[1])
            push!(intervals, (t0, t1, a, b, c, d))
            i += 5; continue
        end
        i += 1
    end

    isempty(intervals) && error("No intervals found in $filepath")

    t_lo = intervals[1][1]
    t_hi = intervals[end][2]

    function phi_prime(t)
        tc  = clamp(t, t_lo, t_hi)
        idx = length(intervals)
        for k in eachindex(intervals)
            if tc <= intervals[k][2]; idx = k; break; end
        end
        t0, _, a, b, c, d = intervals[idx]
        dt = tc - t0
        exp(a + b*dt + c*dt^2 + d*dt^3)
    end

    return phi_prime, metadata, intervals
end

const ϕ_prime, ϕ_meta, ϕ_intervals = parse_polynome(poly_file)

tmin   = ϕ_intervals[1][1]
tmax   = ϕ_intervals[end][2]
t_plot = range(tmin, tmax, length=5000)


## ===== Optimisation helpers =====

# ── decode_params ────────────────────────────────────────────────────────────
# Converts the raw optimiser vector θ_opt to the natural-space parameter
# vector p passed to ODE_MODEL! and PHI_equiv.
#
# Two modes, selected by what the model exposes :
#
#   LOG_INDICES defined  →  new behaviour: exp() applied only to those indices.
#                           All other parameters are passed through unchanged.
#                           (Used by SRVi_phi : α in log, θ_k in linear log-space)
#
#   LOG_INDICES absent   →  legacy behaviour: if LOG_FIT == true, exp() on all;
#                           if LOG_FIT == false, identity.

function make_bounds()
    lower = [fp.lower for fp in ModelDef.FITTED_PARAMS]
    upper = [fp.upper for fp in ModelDef.FITTED_PARAMS]
    return lower, upper
end

# Convert optimiser vector θ → natural-space parameter vector p
function decode_params(θ)
    if isdefined(ModelDef, :LOG_INDICES)
        p = copy(θ)
        for i in ModelDef.LOG_INDICES
            p[i] = exp(θ[i])
        end
        return p
    else
        return ModelDef.LOG_FIT ? exp.(θ) : θ
    end
end

const lower_bounds, upper_bounds = make_bounds()

const isoutofdomain = (u, p, t) -> any(x -> x < 0, u)


## ===== ODE solver =====
function solve_cycle(p, u0, t0, t1, t_save)
    p_ode = (α = p[1], φ = ModelDef.build_phi_func(p))
    prob = ODEProblem(ModelDef.ODE_MODEL!, u0, (t0, t1), p_ode)
    sol  = solve(prob, Rodas5(),
                 reltol=1e-6, abstol=1e-6,
                 saveat=t_save,
                 isoutofdomain=isoutofdomain)
    return sol
end


## ===== State index helpers =====
const _HOST_INDICES  = unique([i for sl in ModelDef.STATE_LABELS
                            for i in sl.indices
                            if sl.state in (:S, :R, :H_total)])
const _VIRUS_INDICES = unique([i for sl in ModelDef.STATE_LABELS
                            for i in sl.indices
                            if sl.state in (:Vi, :Vd, :V_total)])
const _S_INDICES     = unique([i for sl in ModelDef.STATE_LABELS
                            for i in sl.indices
                            if sl.state == :S])
const _VI_INDICES    = unique([i for sl in ModelDef.STATE_LABELS
                            for i in sl.indices
                            if sl.state == :Vi])

_sum_indices(u, idxs) = sum(u[i] for i in idxs)

host_sum(u)  = _sum_indices(u, _HOST_INDICES)
virus_sum(u) = _sum_indices(u, _VIRUS_INDICES)
s_val(u)     = isempty(_S_INDICES)  ? host_sum(u)  : _sum_indices(u, _S_INDICES)
vi_val(u)    = isempty(_VI_INDICES) ? virus_sum(u) : _sum_indices(u, _VI_INDICES)


# ── phi_equiv_at ─────────────────────────────────────────────────────────────
# Unified call that works for both 2-argument (SRVi) and 3-argument (SRVi_phi)
# signatures of ModelDef.PHI_equiv.
function phi_equiv_at(Y_t, p, t)
    m = methods(ModelDef.PHI_equiv)
    # Check the maximum number of positional arguments across all methods
    if any(method.nargs >= 4 for method in m)   # nargs counts 'self' too
        return ModelDef.PHI_equiv(Y_t, p, t)
    else
        return ModelDef.PHI_equiv(Y_t, p)
    end
end


## ===== Objective function =====
function objective(θ; λ_val=λ)
    p = decode_params(θ)

    total_err = 0.0
    prop_S  = 1.0
    prop_Vi = 1.0

    for cycle in 1:cycles_fit
        H0 = mean(raw_data[rep][cycle].H[1] for rep in replicates)
        V0 = mean(raw_data[rep][cycle].V[1] for rep in replicates)
        u0 = ModelDef.INITIAL_CONDITION(H0, V0, prop_S, prop_Vi)

        times = Float64[]
        for rep in replicates
            append!(times, raw_data[rep][cycle].tH)
            append!(times, raw_data[rep][cycle].tV)
        end
        t_data     = sort(unique(times))
        t0_c, t1_c = t_data[1], t_data[end]

        t_phi_grid = range(t0_c, t1_c, length=100)
        t_all      = sort(unique(vcat(t_data, collect(t_phi_grid))))

        sol = solve_cycle(p, u0, t0_c, t1_c, t_all)

        if sol.retcode != SciMLBase.ReturnCode.Success || any(x -> x < 0, reduce(vcat, sol.u))
            return 1e15
        end

        u_end   = sol.u[end]
        H_end   = host_sum(u_end)
        V_end   = virus_sum(u_end)
        prop_S  = s_val(u_end)  / H_end
        prop_Vi = vi_val(u_end) / V_end

        # --- Data residuals ---
        for rep in replicates
            data = raw_data[rep][cycle]

            H_pred = [max(host_sum(sol(t)),  1e-12) for t in data.tH]
            V_pred = [max(virus_sum(sol(t)), 1e-12) for t in data.tV]

            total_err += sum((log.(H_pred) .- log.(data.H)).^2) / length(data.H)
            total_err += sum((log.(V_pred) .- log.(data.V)).^2) / length(data.V)
        end

        # --- φ coherence term ---
        if λ_val > 0
            phi_apparent = [phi_equiv_at(sol(t), p, t)  for t in t_phi_grid]
            phi_ref      = ϕ_prime.(t_phi_grid)
            err_phi      = sum((log.(max.(phi_apparent, 1e-300)) .- log.(phi_ref)).^2) / length(t_phi_grid)
            total_err   += λ_val * err_phi
        end
    end

    return total_err
end


## ===== Optimisation =====
function run_DE(seed)
    Random.seed!(seed)
    res = bboptimize(
        objective;
        SearchRange    = collect(zip(lower_bounds, upper_bounds)),
        NumDimensions  = length(lower_bounds),
        Method         = :xnes,
        PopulationSize = 1000,
        MaxSteps       = 10_000,
        TraceMode      = :silent,
    )
    return (fitness=best_fitness(res), θ=best_candidate(res))
end

function main()
    results = Vector{NamedTuple}(undef, n_runs)
    Threads.@threads for i in 1:n_runs
        println("  Thread $(Threads.threadid()) — run $i (seed=$(1000+i))")
        results[i] = run_DE(1000 + i)
    end
    return results
end

results = Base.invokelatest(main)

best_idx    = argmin(r.fitness for r in results)
best_result = results[best_idx]
θbest       = best_result.θ
p_best      = decode_params(θbest)

println("\n===== Best result ($(ModelDef.MODEL_NAME)) =====")
println("  Loss = ", best_result.fitness)
for (i, fp) in enumerate(ModelDef.FITTED_PARAMS)
    println("  $(fp.name) = $(p_best[i])  ($(fp.description))")
end
println("  All run fitnesses: ", [r.fitness for r in results])


## ===== Simulation & Plot =====

ytick_vals1   = [10.0^i for i in -2:2:8]
ytick_labels1 = [L"10^{%$i}" for i in -2:2:8]
ytick_vals2   = [10.0^i for i in 2:2:10]
ytick_labels2 = [L"10^{%$i}" for i in 2:2:10]

pl_fit = plot(
    layout        = (3, 1),
    size          = (1800, 1600),
    left_margin   = 15mm,
    right_margin  = 10mm,
    top_margin    = 5mm,
    bottom_margin = 10mm,
    grid          = true,
    yscale        = :log10,
    xlims         = (0, 67),
    ytickfontsize  = 22,
    legendfontsize = 15,
    guidefontsize  = 20,
    xtickfontsize  = 20,
    titlefontsize  = 20,
    xlabel         = "Time (days)",
    legend         = :bottomright,
    plot_title     = "Model $(ModelDef.MODEL_NAME) - λ = $(λ)",
    plot_titlefontsize = 25,
)

for rep in replicates
    for cyc_idx in 1:cycles_sim
        data = raw_data[rep][cyc_idx]
        lbl  = cyc_idx == 1 ? "Replicate $rep" : ""
        scatter!(pl_fit[1], data.tH, data.H;
            label=lbl, color=replicate_colors[rep], alpha=0.7,
            ylabel="Host\n(cell/mL)",
            ylims=(1e-2, 3e8), yticks=(ytick_vals1, ytick_labels1),
            markershape=:circle, markersize=8)
        scatter!(pl_fit[2], data.tV, data.V;
            label=lbl, color=replicate_colors[rep], alpha=0.7,
            ylabel="Virus\n(part/mL)",
            ylims=(1e2, 1e10), yticks=(ytick_vals2, ytick_labels2),
            markershape=:circle, markersize=8, legend=(0.13, 0.45))
    end
end

global prop_S_plot  = 1.0
global prop_Vi_plot = 1.0

for cycle in 1:cycles_sim
    H0 = mean(raw_data[rep][cycle].H[1] for rep in replicates)
    V0 = mean(raw_data[rep][cycle].V[1] for rep in replicates)
    u0 = ModelDef.INITIAL_CONDITION(H0, V0, prop_S_plot, prop_Vi_plot)

    t0 = minimum(raw_data[rep][cycle].tH[1]   for rep in replicates)
    t1 = maximum(max(raw_data[rep][cycle].tH[end],
                     raw_data[rep][cycle].tV[end]) for rep in replicates)

    p_plot = (
        α = p_best[1],
        φ = ModelDef.build_phi_func(p_best)
    )
    sol = solve(ODEProblem(ModelDef.ODE_MODEL!, u0, (t0, t1), p_plot),
                Rodas5(), reltol=1e-6, abstol=1e-6,
                isoutofdomain=isoutofdomain)

    u_end        = sol.u[end]
    H_end        = host_sum(u_end)
    V_end        = virus_sum(u_end)
    global prop_S_plot  = H_end > 0 ? clamp(s_val(u_end)  / H_end, 0.0, 1.0) : 1.0
    global prop_Vi_plot = V_end > 0 ? clamp(vi_val(u_end) / V_end, 0.0, 1.0) : 1.0

    lbl(s) = cycle == 1 ? s : ""

    sl_vals(sl) = length(sl.indices) == 1 ? sol[sl.indices[1], :] :
                                             sum(sol[i, :] for i in sl.indices)

    for sl in ModelDef.STATE_LABELS
        sl.state in (:H_total, :S, :R) || continue
        plot!(pl_fit[1], sol.t, sl_vals(sl);
              label=lbl(sl.label), color=sl.color, lw=sl.lw, ls=sl.ls)
    end

    for sl in ModelDef.STATE_LABELS
        sl.state in (:V_total, :Vi, :Vd) || continue
        plot!(pl_fit[2], sol.t, sl_vals(sl);
              label=lbl(sl.label), color=sl.color, lw=sl.lw, ls=sl.ls)
    end

    # Panel 3 : φ_equiv
    phi_app = [phi_equiv_at(sol(t), p_best, t) for t in sol.t]
    plot!(pl_fit[3], sol.t, phi_app;
          label=lbl("φ_equiv"), color=:red, lw=3,
          ylabel="φ_equiv and φ'\n(mL/(cell.day))",
          ylims=(1e-14, 1e-6),
          title="φ_equiv vs φ'")
end

plot!(pl_fit[3], t_plot, ϕ_prime.(t_plot);
      label="φ'", color=:blue, lw=2, ls=:dot,
      xlabel="Time (days)")

display(pl_fit)

mkpath(joinpath(@__DIR__, "output/$(ModelDef.MODEL_NAME)/$(cycles_fit)"))
fig_path = joinpath(@__DIR__, "output/$(ModelDef.MODEL_NAME)/$(cycles_fit)/lambda$(λ).png")
savefig(pl_fit, fig_path)
println("\nFigure saved to $fig_path")


## ===== Log =====
log_path = joinpath(@__DIR__, "output/$(ModelDef.MODEL_NAME)/$(cycles_fit)/lambda$(λ)_log.txt")
open(log_path, "w") do io
    println(io, "===== Run log =====")
    println(io, "Date        : ", Dates.now())
    println(io, "Model       : ", ModelDef.MODEL_NAME)
    println(io, "Model file  : ", MODEL_FILE)
    println(io, "λ           : ", λ)
    println(io, "n_runs      : ", n_runs)
    println(io, "cycles_fit  : ", cycles_fit)
    println(io, "cycles_sim  : ", cycles_sim)
    println(io, "replicates  : ", replicates)
    println(io, "ϕ_combi     : ", ϕ_combi)
    println(io)
    println(io, "===== Fixed parameters =====")
    for (k, v) in pairs(ModelDef.FIXED_PARAMS)
        println(io, "  $k = $v")
    end
    if isdefined(ModelDef, :KNOTS)
        println(io, "  KNOTS = ", ModelDef.KNOTS)
    end
    println(io)
    println(io, "===== Search bounds =====")
    println(io, "  LOG_FIT     : ", ModelDef.LOG_FIT)
    if isdefined(ModelDef, :LOG_INDICES)
        println(io, "  LOG_INDICES : ", ModelDef.LOG_INDICES)
    end
    for fp in ModelDef.FITTED_PARAMS
        println(io, "  $(fp.name)  ∈ [$(fp.lower), $(fp.upper)]  — $(fp.description)")
    end
    println(io)
    println(io, "===== All runs =====")
    for (i, res) in enumerate(results)
        p_i = decode_params(res.θ)
        param_str = join(["$(ModelDef.FITTED_PARAMS[j].name)=$(p_i[j])" for j in eachindex(ModelDef.FITTED_PARAMS)], "  ")
        println(io, "  Run $i (seed=$(1000+i)) : loss=$(res.fitness)  $param_str")
    end
    println(io)
    println(io, "===== Best result =====")
    println(io, "  Loss = ", best_result.fitness)
    for (i, fp) in enumerate(ModelDef.FITTED_PARAMS)
        println(io, "  $(fp.name) = $(p_best[i])  ($(fp.description))")
    end
    println(io)
    println(io, "===== Polynomial =====")
    println(io, "  File        : ", poly_file)
    println(io, "  n_intervals : ", length(ϕ_intervals))
    haskey(ϕ_meta, "best_fitness") && println(io, "  best_fitness : ", ϕ_meta["best_fitness"])
end
println("Log saved to $log_path")


## ===== S/R ratio diagnostic =====
println("\n===== S/R ratio along cycles =====")
global prop_S_diag  = 1.0
global prop_Vi_diag = 1.0

for cycle in 1:cycles_sim
    H0 = mean(raw_data[rep][cycle].H[1] for rep in replicates)
    V0 = mean(raw_data[rep][cycle].V[1] for rep in replicates)
    u0 = ModelDef.INITIAL_CONDITION(H0, V0, prop_S_diag, prop_Vi_diag)

    t0 = minimum(raw_data[rep][cycle].tH[1]   for rep in replicates)
    t1 = maximum(max(raw_data[rep][cycle].tH[end],
                     raw_data[rep][cycle].tV[end]) for rep in replicates)
    p_plot = (
        α = p_best[1],
        φ = ModelDef.build_phi_func(p_best)
    )
    sol = solve(ODEProblem(ModelDef.ODE_MODEL!, u0, (t0, t1), p_plot),
                Rodas5(), reltol=1e-6, abstol=1e-6,
                isoutofdomain=isoutofdomain)

    println("  Cycle $cycle (t = $t0 → $t1 days)")
    println("    u0 : S=$(u0[1])  R=$(u0[2])  prop_S=$(prop_S_diag)")
    for t in range(t0, t1, length=6)
        u  = sol(t)
        S  = u[1];  R = u[2]
        ratio = R > 1e-30 ? S / R : Inf
        println("    t=$(round(t, digits=1)) d :  S=$(round(S, sigdigits=3))  R=$(round(R, sigdigits=3))  S/R=$(round(ratio, sigdigits=3))  S/H=$(round(S/(S+R), sigdigits=3))")
    end

    u_end        = sol.u[end]
    H_end        = host_sum(u_end)
    V_end        = virus_sum(u_end)
    global prop_S_diag  = H_end > 0 ? clamp(s_val(u_end)  / H_end, 0.0, 1.0) : 1.0
    global prop_Vi_diag = V_end > 0 ? clamp(vi_val(u_end) / V_end, 0.0, 1.0) : 1.0
    println("    → end of cycle : prop_S=$(round(prop_S_diag, sigdigits=3))  prop_Vi=$(round(prop_Vi_diag, sigdigits=3))")
end