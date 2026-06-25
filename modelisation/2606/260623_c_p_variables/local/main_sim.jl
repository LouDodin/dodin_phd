## ===== Packages =====
using DifferentialEquations
using OrdinaryDiffEqRosenbrock
using CSV
using DataFrames
using Statistics
using Plots
using Measures
using LaTeXStrings
using SciMLBase
using Colors
using DataInterpolations


## ===== Choices =====
const MODEL_FILE = joinpath(@__DIR__, "models/SR_RS.jl")
const replicates = ["A", "B", "C"]
const cycles_sim = 5

const LOG_FILE = joinpath(@__DIR__, "output/SR_RS/1_1_1_1_1/SR_RS_spline_log.txt")


## ===== Utils =====
include(MODEL_FILE)
using .ModelDef

const isoutofdomain = (u, p, t) -> any(x -> x < 0, u)

replicate_colors = Dict(
    "A" => RGB(0.6, 0.8, 1.0),
    "B" => RGB(31/255, 119/255, 180/255),
    "C" => RGB(0.0, 0.3, 0.7)
)
model_color = RGB(255/255, 127/255, 14/255)
cycle_color = RGB(31/255, 119/255, 180/255)


## ===== Input: H & V data =====

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

const t_global_min = minimum(raw_data[rep][1].tH[1] for rep in replicates)
const t_global_max = maximum(max(raw_data[rep][cycles_sim].tH[end],
                                 raw_data[rep][cycles_sim].tV[end]) for rep in replicates)


## ===== Log parser =====

function parse_spline_log(filepath::String)
    isfile(filepath) || error("Log file not found: $filepath")
    lines = readlines(filepath)

    ϕ_val      = nothing
    c_vals     = Float64[]
    p_vals     = Float64[]
    n_interior = Int[]
    fitness    = nothing
    aic_val    = nothing
    aicc_val   = nothing
    bic_val    = nothing

    for line in lines
        m = match(r"^\s*ϕ\s*=\s*([0-9eE+\-.]+)", line)
        if m !== nothing; ϕ_val = parse(Float64, m.captures[1]); end

        m = match(r"c at knots\s*:\s*\[(.+)\]", line)
        if m !== nothing; c_vals = parse.(Float64, split(m.captures[1], ',')); end

        m = match(r"p at knots\s*:\s*\[(.+)\]", line)
        if m !== nothing; p_vals = parse.(Float64, split(m.captures[1], ',')); end

        m = match(r"n_interior\s*:\s*\[(.+)\]", line)
        if m !== nothing; n_interior = parse.(Int, split(m.captures[1], ',')); end

        m = match(r"Loss\s*=\s*([0-9eE+\-.]+)", line)
        if m !== nothing; fitness = parse(Float64, m.captures[1]); end

        m = match(r"^\s*AIC\s*=\s*([0-9eE+\-.]+)",  line); m !== nothing && (aic_val  = parse(Float64, m.captures[1]))
        m = match(r"^\s*AICc\s*=\s*([0-9eE+\-.]+)", line); m !== nothing && (aicc_val = parse(Float64, m.captures[1]))
        m = match(r"^\s*BIC\s*=\s*([0-9eE+\-.]+)",  line); m !== nothing && (bic_val  = parse(Float64, m.captures[1]))
    end

    ϕ_val   === nothing && error("Could not parse ϕ")
    isempty(c_vals)     && error("Could not parse c knots")
    isempty(p_vals)     && error("Could not parse p knots")
    isempty(n_interior) && error("Could not parse n_interior")

    println("Parsed from log:")
    println("  ϕ         = $ϕ_val")
    println("  N_KNOTS   = $(length(c_vals))")
    println("  fitness   = $fitness")
    println("  AIC=$aic_val  AICc=$aicc_val  BIC=$bic_val")

    return (ϕ=ϕ_val, c_vals=c_vals, p_vals=p_vals, n_interior=n_interior,
            fitness=fitness, AIC=aic_val, AICc=aicc_val, BIC=bic_val)
end

parsed = parse_spline_log(LOG_FILE)

const ϕ_best  = parsed.ϕ
const N_KNOTS = length(parsed.c_vals)


## ===== Spline knot reconstruction =====

function build_spline_knots(raw_data, replicates, cycles_fit, n_interior)
    knots = Float64[]
    for cyc in 1:cycles_fit
        t0 = minimum(raw_data[rep][cyc].tH[1] for rep in replicates)
        t1 = maximum(max(raw_data[rep][cyc].tH[end],
                         raw_data[rep][cyc].tV[end]) for rep in replicates)
        push!(knots, t0)
        for j in 1:n_interior[cyc]
            push!(knots, t0 + j * (t1 - t0) / (n_interior[cyc] + 1))
        end
    end
    cyc = length(n_interior)
    t1_last = maximum(max(raw_data[rep][cyc].tH[end],
                          raw_data[rep][cyc].tV[end]) for rep in replicates)
    push!(knots, t1_last)
    return knots
end

const cycles_fit   = length(parsed.n_interior)
const spline_knots = build_spline_knots(raw_data, replicates, cycles_fit, parsed.n_interior)

length(spline_knots) == N_KNOTS ||
    @warn "Knot count mismatch: rebuilt $(length(spline_knots)) vs log $N_KNOTS"

println("Spline knots reconstructed ($(length(spline_knots))): ", round.(spline_knots, digits=2))


## ===== Spline reconstruction =====

sigmoid(x) = 1.0 / (1.0 + exp(-x))

function build_spline_func(knots::Vector{Float64}, values::Vector{Float64})
    itp = CubicSpline(values, knots)
    return t -> itp(clamp(t, knots[1], knots[end]))
end

log_c_func   = build_spline_func(spline_knots, log.(parsed.c_vals))
logit_p_func = build_spline_func(spline_knots, log.(parsed.p_vals ./ (1.0 .- parsed.p_vals)))

c_func_best = t -> exp(log_c_func(t))
p_func_best = t -> sigmoid(logit_p_func(t))


## ===== ϕ_ref =====

const poly_file = "modelisation/2606/260604/genotoul/output/nint_3_2_2_3_2/polynomial.txt"

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

    t_lo = intervals[1][1];  t_hi = intervals[end][2]

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

    return phi_ref_func, knots_ref
end

const ϕ_ref_func, knots_ref = parse_phi_ref(poly_file)

const N_GRID     = 10_000
const t_grid     = collect(range(t_global_min, t_global_max; length=N_GRID))
const ϕ_ref_grid = ϕ_ref_func.(t_grid)


## ===== Plot =====

ytick_vals1   = [10.0^i for i in 2:9]
ytick_labels1 = [L"10^{%$i}" for i in 2:9]
ytick_vals2   = [10.0^i for i in 5:10]
ytick_labels2 = [L"10^{%$i}" for i in 5:10]
ytick_vals3   = [10.0^i for i in -14:-7]
ytick_labels3 = [L"10^{%$i}" for i in -14:-7]

pl_sim = plot(
    layout             = (3, 1),
    size               = (1800, 2000),
    left_margin        = 15mm,
    right_margin       = 10mm,
    top_margin         = 35mm,
    bottom_margin      = 5mm,
    grid               = true,
    yscale             = :log10,
    xlims              = (t_global_min, t_global_max),
    ytickfontsize      = 22,
    legendfontsize     = 15,
    guidefontsize      = 20,
    xtickfontsize      = 20,
    xlabel             = "Time [days]",
    legend             = :bottomright,
    plot_title         = "$(ModelDef.MODEL_NAME)\nFitness = $(round(something(parsed.fitness, NaN), sigdigits=4))\nAIC = $(round(something(parsed.AIC, NaN), sigdigits=4)) - AICc = $(round(something(parsed.AICc, NaN), sigdigits=4)) - BIC = $(round(something(parsed.BIC, NaN), sigdigits=4))",
    plot_titlefontsize = 25,
)

# ── Data ─────────────────────────────────────────────────────────────────
for rep in replicates
    for cyc_idx in 1:cycles_sim
        data = raw_data[rep][cyc_idx]
        lbl  = cyc_idx == 1 ? "Replicate $rep" : ""
        scatter!(pl_sim[1], data.tH, data.H;
            label=lbl, color=replicate_colors[rep], alpha=0.7,
            ylabel="Host\n[cell/mL]",
            ylims=(1e2, 1e9), yticks=(ytick_vals1, ytick_labels1),
            markershape=:circle, markersize=10,
            legend=(0.15, 1))
        scatter!(pl_sim[2], data.tV, data.V;
            label=lbl, color=replicate_colors[rep], alpha=0.7,
            ylabel="Virus\n[part/mL]",
            ylims=(1e5, 1e10), yticks=(ytick_vals2, ytick_labels2),
            markershape=:circle, markersize=10,
            legend=(0.15, 0.3))
    end
end

plot!(pl_sim[3], t_grid, ϕ_ref_grid;
    label="ϕ_ref", color=replicate_colors["B"], lw=6, alpha=0.7,
    ylabel="ϕ\n[mL/(part.day)]",
    ylims=(1e-14, 1e-7), yticks=(ytick_vals3, ytick_labels3))

scatter!(pl_sim[3], knots_ref, ϕ_ref_func.(knots_ref);
    color=replicate_colors["C"], markersize=10, markershape=:diamond,
    label="Knots", legend=:bottomright)

# ── Model trajectories ────────────────────────────────────────────────────
global prop_S_plot = 1.0

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

    u_end = sol.u[end];  H_end = u_end[1] + u_end[2]
    global prop_S_plot = H_end > 0 ? clamp(u_end[1] / H_end, 0.0, 1.0) : 1.0

    S = [u[1] for u in sol.u]
    R = [u[2] for u in sol.u]
    V = [u[3] for u in sol.u]
    H = S .+ R
    ϕ_equiv = [ModelDef.PHI_EQUIV(sol(t), (c_func_best, p_func_best, ϕ_best)) for t in sol.t]

    lbl1 = cycle == 1
    plot!(pl_sim[1], sol.t, H; label=lbl1 ? "H model" : "", color=model_color,  lw=6, alpha=0.7)
    plot!(pl_sim[1], sol.t, S; label=lbl1 ? "S model" : "", color=:red,          lw=4, ls=:dash, alpha=0.7)
    plot!(pl_sim[1], sol.t, R; label=lbl1 ? "R model" : "", color=:green,        lw=4, ls=:dash, alpha=0.7)
    plot!(pl_sim[2], sol.t, V; label=lbl1 ? "V model" : "", color=model_color,   lw=6, alpha=0.7)
    plot!(pl_sim[3], sol.t, ϕ_equiv; label=lbl1 ? "ϕ_equiv model" : "", color=model_color, lw=6, alpha=0.7)
end

# ── Dilution lines ────────────────────────────────────────────────────────
cycle_tbounds = [(
    minimum(raw_data[rep][cyc].tH[1] for rep in replicates),
    maximum(max(raw_data[rep][cyc].tH[end], raw_data[rep][cyc].tV[end]) for rep in replicates)
) for cyc in 1:cycles_sim]

dilution_times = [cycle_tbounds[c][2] for c in 1:cycles_sim-1]

for panel in 1:3
    for td in dilution_times
        vline!(pl_sim[panel], [td]; color=cycle_color, lw=2, ls=:dot, label=nothing)
    end
end

display(pl_sim)

fig_path = joinpath(dirname(LOG_FILE), "$(ModelDef.MODEL_NAME)_sim_plot.png")
savefig(pl_sim, fig_path)
println("\nFigure saved to $fig_path")