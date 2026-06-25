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
const MODEL_FILE  = joinpath(@__DIR__, "models/SR_RS.jl")
const replicates  = ["A", "B", "C"]
const cycles_fit  = 5
const cycles_sim  = 5
const n_runs      = 1
const n_interior  = [2,2,2,2,2]   # interior knots per cycle for c and p splines


## ===== Utils =====
include(MODEL_FILE)
using .ModelDef

output_dir = joinpath(@__DIR__, "output/$(ModelDef.MODEL_NAME)/$(n_interior[1])_$(n_interior[2])_$(n_interior[3])_$(n_interior[4])_$(n_interior[5])")
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
    ϕ       = exp(θ[1])
    log_c   = θ[2 : N_KNOTS+1]
    logit_p = θ[N_KNOTS+2 : 2*N_KNOTS+1]

    # Splines construites en espace transformé (log / logit)
    log_c_func   = build_spline_func(spline_knots, log_c)
    logit_p_func = build_spline_func(spline_knots, logit_p)

    # Fonctions finales : transformation appliquée APRÈS interpolation
    c_func = t -> exp(log_c_func(t))
    p_func = t -> sigmoid(logit_p_func(t))

    # Valeurs aux nœuds pour le plot
    c_vals = exp.(log_c)
    p_vals = sigmoid.(logit_p)

    return ϕ, c_func, p_func, c_vals, p_vals
end


## ===== Search bounds =====
#
# From ModelDef.FITTED_PARAMS:
#   :ϕ  → log-space bounds for θ[1]
#   :c  → log-space bounds for θ[2..N_KNOTS+1]
#   :p  → logit-space bounds for θ[N_KNOTS+2..2*N_KNOTS+1]

fp_dict = Dict(fp.name => fp for fp in ModelDef.FITTED_PARAMS)

log_ϕ_lo  = log(fp_dict[:ϕ].lower);   log_ϕ_hi  = log(fp_dict[:ϕ].upper)
log_c_lo  = log(fp_dict[:c].lower);   log_c_hi  = log(fp_dict[:c].upper)
logit_p_lo = log(fp_dict[:p].lower / (1 - fp_dict[:p].lower))
logit_p_hi = log(fp_dict[:p].upper / (1 - fp_dict[:p].upper))

search_range = vcat(
    [(log_ϕ_lo, log_ϕ_hi)],
    fill((log_c_lo, log_c_hi), N_KNOTS),
    fill((logit_p_lo, logit_p_hi), N_KNOTS)
)

const N_PARAMS = 1 + 2 * N_KNOTS
println("Total parameters: $(N_PARAMS)  (1 ϕ + $(N_KNOTS) c-knots + $(N_KNOTS) p-knots)")


## ===== ODE solver =====

function solve_cycle(c_func, p_func, ϕ, u0, t0, t1, t_save)
    params = (c_func, p_func, ϕ)
    prob   = ODEProblem(ModelDef.ODE_MODEL!, u0, (t0, t1), params)
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
    ϕ, c_func, p_func, _, _ = decode_theta(θ)

    total_err = 0.0

    for rep in replicates
        prop_S = 1.0

        for cycle in 1:cycles_fit
            data = raw_data[rep][cycle]
            H0   = data.H[1];  V0 = data.V[1]
            u0   = ModelDef.INITIAL_CONDITION(H0, V0, prop_S)

            t_data = sort(unique(vcat(data.tH, data.tV)))
            t0_data, t1_data = t_data[1], t_data[end]

            sol = solve_cycle(c_func, p_func, ϕ, u0, t0_data, t1_data, t_data)

            if sol.retcode != SciMLBase.ReturnCode.Success ||
               any(x -> x < 0, reduce(vcat, sol.u))
                return 1e15
            end

            H_pred = [sol(t)[1] + sol(t)[2] for t in data.tH]
            V_pred = [sol(t)[3]              for t in data.tV]

            total_err += sum((log.(H_pred) .- log.(data.H)).^2) / length(data.H)
            total_err += sum((log.(V_pred) .- log.(data.V)).^2) / length(data.V)

            u_end  = sol.u[end]
            H_end  = u_end[1] + u_end[2]
            prop_S = H_end > 0 ? clamp(u_end[1] / H_end, 0.0, 1.0) : 1.0
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
        NumDimensions        = N_PARAMS,
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
best_idx  = argmin(r.fitness for r in results)
best_result = results[best_idx]
θbest     = best_result.θ

ϕ_best, c_func_best, p_func_best, c_knot_vals, p_knot_vals = decode_theta(θbest)

println("\n===== Best result ($(ModelDef.MODEL_NAME)) =====")
println("  Fitness = ", best_result.fitness)
println("  ϕ       = ", ϕ_best)
println("  c at knots: ", round.(c_knot_vals, sigdigits=3))
println("  p at knots: ", round.(p_knot_vals, sigdigits=3))


## ===== Information Criteria =====

n_obs = sum(
    length(raw_data[rep][cyc].H) + length(raw_data[rep][cyc].V)
    for rep in replicates
    for cyc in 1:cycles_fit
)

function raw_rss(θ)
    ϕ, c_func, p_func, _, _ = decode_theta(θ)
    rss = 0.0

    for rep in replicates
        prop_S = 1.0
        for cycle in 1:cycles_fit
            data = raw_data[rep][cycle]
            H0   = data.H[1];  V0 = data.V[1]
            u0   = ModelDef.INITIAL_CONDITION(H0, V0, prop_S)

            t_data = sort(unique(vcat(data.tH, data.tV)))
            sol    = solve_cycle(c_func, p_func, ϕ, u0, t_data[1], t_data[end], t_data)

            (sol.retcode != SciMLBase.ReturnCode.Success ||
             any(x -> x < 0, reduce(vcat, sol.u))) && return Inf

            H_pred = [sol(t)[1] + sol(t)[2] for t in data.tH]
            V_pred = [sol(t)[3]              for t in data.tV]

            rss += sum((log.(H_pred) .- log.(data.H)).^2)
            rss += sum((log.(V_pred) .- log.(data.V)).^2)

            u_end  = sol.u[end]
            H_end  = u_end[1] + u_end[2]
            prop_S = H_end > 0 ? clamp(u_end[1] / H_end, 0.0, 1.0) : 1.0
        end
    end
    return rss
end

rss_best = raw_rss(θbest)
k        = N_PARAMS
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


## ===== ϕ_ref (reference polynomial for comparison) =====

const phi_ref_combi = "3_2_2_3_2"
const phi_ref_file  = joinpath(@__DIR__,
    "../../260604/genotoul/output/nint_$(phi_ref_combi)/polynomial.txt")

function parse_phi_ref(filepath::String)
    isfile(filepath) || error("ϕ_ref file not found: $filepath")
    lines = readlines(filepath)

    function extract_float(s::AbstractString)::Float64
        m = match(r"([+-]?\s*[0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)", strip(s))
        m === nothing && error("Cannot extract float from: \"$s\"")
        parse(Float64, replace(m.captures[1], " " => ""))
    end

    intervals = NTuple{6,Float64}[]
    i = 1
    while i <= length(lines)
        m_iv = match(r"Interval\s+\[([0-9eE+\-.]+),\s*([0-9eE+\-.]+)\]\s+days:", strip(lines[i]))
        if m_iv !== nothing
            if i + 4 > length(lines); i += 1; continue; end
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

    phi_ref_func = function(t)
        tc  = clamp(t, t_lo, t_hi)
        idx = length(intervals)
        for k in eachindex(intervals)
            if tc <= intervals[k][2]; idx = k; break; end
        end
        t0, _, a, b, c_coef, d = intervals[idx]
        dt = tc - t0
        exp(a + b*dt + c_coef*dt^2 + d*dt^3)
    end

    knots_ref = Float64[]
    for iv in intervals; push!(knots_ref, iv[1]); end
    push!(knots_ref, intervals[end][2])
    unique!(sort!(knots_ref))

    return phi_ref_func, intervals, knots_ref
end

const ϕ_ref_func, phi_ref_intervals, knots_ref = parse_phi_ref(phi_ref_file)
println("ϕ_ref parsed: $(length(phi_ref_intervals)) intervals")

const N_GRID  = 10_000
const t_global_min = minimum(raw_data[rep][1].tH[1]                                  for rep in replicates)
const t_global_max = maximum(max(raw_data[rep][cycles_sim].tH[end],
                                 raw_data[rep][cycles_sim].tV[end]) for rep in replicates)
const t_grid  = collect(range(t_global_min, t_global_max; length=N_GRID))
const ϕ_ref_grid = ϕ_ref_func.(t_grid)


## ===== Simulation & Plot =====

ytick_vals_H   = [10.0^i for i in 2:9]
ytick_labels_H = [L"10^{%$i}" for i in 2:9]
ytick_vals_V   = [10.0^i for i in 5:10]
ytick_labels_V = [L"10^{%$i}" for i in 5:10]
ytick_vals_phi   = [10.0^i for i in -20:-7]
ytick_labels_phi = [L"10^{%$i}" for i in -20:-7]

ytick_vals_c   = [10.0^i for i in -9:0]
ytick_labels_c = [L"10^{%$i}" for i in -9:0]

ytick_vals_p   = 0:0.1:1
ytick_labels_p = string.(ytick_vals_p)

ytick_vals_alpha   = [10.0^i for i in -15:0]
ytick_labels_alpha = [L"10^{%$i}" for i in -15:0]

ytick_vals_gamma   = [10.0^i for i in -15:0]
ytick_labels_gamma = [L"10^{%$i}" for i in -15:0]


title_str = "$(ModelDef.MODEL_NAME)  |  ϕ=$(round(ϕ_best, sigdigits=3))\n" *
            "Fitness=$(round(best_result.fitness, sigdigits=4))  " *
            "AIC=$(round(AIC, sigdigits=4))  AICc=$(round(AICc, sigdigits=4))  BIC=$(round(BIC, sigdigits=4))"

common_kw = (
    left_margin   = 15mm,
    right_margin  = 10mm,
    bottom_margin = 5mm,
    grid          = true,
    xlims         = (t_global_min, t_global_max),
    xtickfontsize = 14,
    ytickfontsize = 14,
    legendfontsize= 11,
    guidefontsize = 15,
    xlabel        = "Time [days]",
    legend        = :bottomright,
)

pl_fit = plot(
    layout            = (4, 2),
    size              = (2500, 2000),
    top_margin        = 10mm,
    plot_title        = title_str,
    plot_titlefontsize= 20,
)

# ── Subplot 1: H, S, R ────────────────────────────────────────────────────
for rep in replicates
    for cyc in 1:cycles_sim
        data = raw_data[rep][cyc]
        scatter!(pl_fit[1], data.tH, data.H;
            color = replicate_colors[rep],
            label = cyc == 1 ? "H Rep $rep" : "",
            ylabel = "Host [cell/mL]",
            yscale = :log10, ylims = (1e2, 1e9),
            yticks = (ytick_vals_H, ytick_labels_H),
            markersize = 7, alpha = 0.7,
            common_kw...,
        )
    end
end

# ── Subplot 2: V ──────────────────────────────────────────────────────────
for rep in replicates
    for cyc in 1:cycles_sim
        data = raw_data[rep][cyc]
        scatter!(pl_fit[3], data.tV, data.V;
            color = replicate_colors[rep],
            label = cyc == 1 ? "V Rep $rep" : "",
            ylabel = "Virus [part/mL]",
            yscale = :log10, ylims = (1e5, 1e10),
            yticks = (ytick_vals_V, ytick_labels_V),
            markersize = 7, alpha = 0.7,
            common_kw...,
        )
    end
end

# ── Subplot 3: ϕ_ref and ϕ_equiv (model) ─────────────────────────────────
plot!(pl_fit[5], t_grid, ϕ_ref_grid;
    label  = "ϕ_ref", color = replicate_colors["B"],
    lw = 3, alpha = 0.8,
    ylabel = "ϕ [mL/(part·day)]",
    yscale = :log10, ylims = (1e-20, 1e-7),
    yticks = (ytick_vals_phi, ytick_labels_phi),
    common_kw...,
)
phi_ref_knot_vals = ϕ_ref_func.(knots_ref)
scatter!(pl_fit[5], knots_ref, phi_ref_knot_vals;
    color = replicate_colors["C"], markersize = 8,
    markershape = :diamond, label = "ϕ_ref knots",
)

# ── Subplots 4–7: c, p, α=c*p, γ=c*(1-p) ─────────────────────────────────
c_grid = c_func_best.(t_grid)
p_grid = p_func_best.(t_grid)
c_ref = 0.0095
p_ref = 0.0023

plot!(pl_fit[2], t_grid, c_grid;
    label = "c", color = model_color, lw = 3,
    ylabel = "c", yscale = :log10, ylims = (1e-9, 1e0), yticks=(ytick_vals_c, ytick_labels_c),  
    common_kw...,
)
scatter!(pl_fit[2], spline_knots, c_knot_vals;
    color = :red, markersize = 8, markershape = :circle, label = "knots",
)
plot!(pl_fit[2], t_grid, fill(c_ref, N_GRID);
    label = "c_ref", color = :black, lw = 3, ls = :dash,
)


plot!(pl_fit[4], t_grid, p_grid;
    label = "p", color = model_color, lw = 3,
    ylabel = "p ∈ (0,1)", ylims = (0-1e-2, 1+1e-2), yticks=(ytick_vals_p, ytick_labels_p),
    common_kw...,
)
scatter!(pl_fit[4], spline_knots, p_knot_vals;
    color = :red, markersize = 8, markershape = :circle, label = "knots",
)
plot!(pl_fit[4], t_grid, fill(p_ref, N_GRID);
    label = "p_ref", color = :black, lw = 3, ls = :dash,
)


plot!(pl_fit[6], t_grid, c_grid .* p_grid;
    label = "α = c·p", color = model_color, lw = 3,
    ylabel = "α [day⁻¹]", yscale = :log10, ylims = (1e-15, 1e0), yticks=(ytick_vals_alpha, ytick_labels_alpha), legend=:topright,
    common_kw...,
)
scatter!(pl_fit[6], spline_knots, c_knot_vals .* p_knot_vals;
    color = :red, markersize = 8, markershape = :circle, label = "knots",
)
plot!(pl_fit[6], t_grid, fill(c_ref * p_ref, N_GRID);
    label = "α_ref", color = :black, lw = 3, ls = :dash,
)


plot!(pl_fit[8], t_grid, c_grid .* (1 .- p_grid);
    label = "γ = c·(1−p)", color = model_color, lw = 3,
    ylabel = "γ [day⁻¹]", yscale = :log10, ylims = (1e-15, 1e0), yticks=(ytick_vals_gamma, ytick_labels_gamma),
    common_kw...,
)
scatter!(pl_fit[8], spline_knots, c_knot_vals .* (1 .- p_knot_vals);
    color = :red, markersize = 8, markershape = :circle, label = "knots",
)
plot!(pl_fit[8], t_grid, fill(c_ref * (1 - p_ref), N_GRID);
    label = "γ_ref", color = :black, lw = 3, ls = :dash,
)



# ── Model trajectories (simulated from best params) ───────────────────────
global prop_S_plot = 1.0
global first_phi   = true

for cycle in 1:cycles_sim
    H0 = mean(raw_data[rep][cycle].H[1] for rep in replicates)
    V0 = mean(raw_data[rep][cycle].V[1] for rep in replicates)
    u0 = ModelDef.INITIAL_CONDITION(H0, V0, prop_S_plot)

    t0 = minimum(raw_data[rep][cycle].tH[1]   for rep in replicates)
    t1 = maximum(max(raw_data[rep][cycle].tH[end],
                     raw_data[rep][cycle].tV[end]) for rep in replicates)

    sol = solve(
        ODEProblem(ModelDef.ODE_MODEL!, u0, (t0, t1), (c_func_best, p_func_best, ϕ_best)),
        Rodas5(), reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain,
    )

    u_end  = sol.u[end]
    H_end  = u_end[1] + u_end[2]
    global prop_S_plot = H_end > 0 ? clamp(u_end[1] / H_end, 0.0, 1.0) : 1.0

    S_traj = [u[1] for u in sol.u]
    R_traj = [u[2] for u in sol.u]
    H_traj = S_traj .+ R_traj
    V_traj = [u[3] for u in sol.u]

    # ϕ_equiv = ϕ * S / H
    phi_equiv_traj = [ModelDef.PHI_EQUIV(sol(t), (c_func_best, p_func_best, ϕ_best))
                      for t in sol.t]

    lbl1 = cycle == 1
    plot!(pl_fit[1], sol.t, H_traj; label=lbl1 ? "H model" : "", color=model_color, lw=4, alpha=0.8)
    plot!(pl_fit[1], sol.t, S_traj; label=lbl1 ? "S model" : "", color=:red,         lw=3, ls=:dash, alpha=0.8)
    plot!(pl_fit[1], sol.t, R_traj; label=lbl1 ? "R model" : "", color=:green,       lw=3, ls=:dash, alpha=0.8)
    plot!(pl_fit[3], sol.t, V_traj; label=lbl1 ? "V model" : "", color=model_color,  lw=4, alpha=0.8)
    plot!(pl_fit[5], sol.t, phi_equiv_traj;
          label = lbl1 ? "ϕ_equiv model" : "",
          color = model_color, lw = 4, alpha = 0.8,
    )
end

# ── Dilution lines on all 7 panels ───────────────────────────────────────
cycle_tbounds = [(
    minimum(raw_data[rep][cyc].tH[1] for rep in replicates),
    maximum(max(raw_data[rep][cyc].tH[end], raw_data[rep][cyc].tV[end]) for rep in replicates)
) for cyc in 1:cycles_sim]

dilution_times = [cycle_tbounds[c][2] for c in 1:cycles_sim-1]

for panel in 1:7
    for td in dilution_times
        vline!(pl_fit[panel], [td]; color=cycle_color, lw=2, ls=:dot, label=nothing)
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
    println(io, "N_PARAMS    : ", N_PARAMS, "  (1 ϕ + $(N_KNOTS) c-knots + $(N_KNOTS) p-knots)")
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
    println(io, "  ϕ    = ", ϕ_best)
    println(io, "  c at knots : ", c_knot_vals)
    println(io, "  p at knots : ", p_knot_vals)
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