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
const ϕ_combi     = "3_2_2_3_2"
const replicates  = ["A", "B", "C"]
const cycles_fit  = 5
const cycles_sim  = 5
const n_runs      = 1

# λ : weight of the φ coherence term
const λ = 0

# Limit display to a specific knot (nothing = no limit, or specify knot index)
const knot_limit = 1#nothing  # Set to knot index (e.g., 3) to limit display


## ===== Colors =====
color_A = RGB(0.6, 0.8, 1.0)
color_B = RGB(31/255, 119/255, 180/255)
color_C = RGB(0.0, 0.3, 0.7)
model_color = RGB(255/255, 127/255, 14/255)
data_color = RGB(31/255, 119/255, 180/255)

replicate_colors = Dict(
    "A" => color_A,
    "B" => color_B,
    "C" => color_C
)


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

# Build interpolation grid for model
tmin   = ϕ_intervals[1][1]
tmax   = ϕ_intervals[end][2]
t_plot = range(tmin, tmax, length=5000)

# Create linear interpolation of phi for ODE model
const N_GRID  = 10_000
const t_grid  = collect(range(tmin, tmax; length=N_GRID))
const phi_grid = ϕ_prime.(t_grid)
const ϕ_interp = LinearInterpolation(phi_grid, t_grid)

# Extract knots from intervals
knots = Float64[]
for interval in ϕ_intervals
    push!(knots, interval[1])
end
if !isempty(ϕ_intervals)
    push!(knots, ϕ_intervals[end][2])
end
unique!(sort!(knots))

# Calculate time limit based on knot_limit
t_limit = if knot_limit === nothing
    tmax
else
    if knot_limit < 1 || knot_limit > length(knots)
        error("knot_limit must be between 1 and $(length(knots))")
    end
    knots[knot_limit]
end
t_plot_limited = range(tmin, min(t_limit, tmax), length=5000)


## ===== Constants =====
r = 0.5592225270686286
K = 7.29695252684594e7
β = 144
δ = 0.02


## ===== Model =====
function model!(dY, Y, p, t)
    S, V = Y[1], Y[2]
    ϕt = ϕ_interp(clamp(t, t_grid[1], t_grid[end]))
    dY[1] = r*S*(1 - S/K) - ϕt*S*V
    dY[2] = β*ϕt*S*V - δ*V
end

isoutofdomain(u, p, t) = any(x -> x < 0, u)


## ===== Simulation & Plot =====

ytick_vals1   = [10.0^i for i in 2:9]
ytick_labels1 = [L"10^{%$i}" for i in 2:9]
ytick_vals2   = [10.0^i for i in 5:10]
ytick_labels2 = [L"10^{%$i}" for i in 5:10]
ytick_vals3   = [10.0^i for i in -13:-8]
ytick_labels3 = [L"10^{%$i}" for i in -13:-8]

pl_fit = plot(
    layout = (3, 1),
    size = (1800, 1500),
    left_margin = 15mm,
    right_margin = 10mm,
    top_margin = 5mm,
    bottom_margin = 10mm,
    grid = true,
    yscale = :log10,
    xlims = (0, 67),
    ytickfontsize = 22,
    legendfontsize = 15,
    guidefontsize = 20,
    xtickfontsize = 20,
    titlefontsize = 20,
    xlabel = "Time [days]",
    legend = :bottomright,
    plot_title = "Model with variable phi (φ: $ϕ_combi) - λ = $(λ)" * (knot_limit === nothing ? "" : " [Knot limit: $knot_limit]"),
    plot_titlefontsize = 25
)

# --- Data & Model for H ---
for rep in replicates
    for cyc_idx in 1:cycles_sim
        data = raw_data[rep][cyc_idx]
        lbl  = cyc_idx == 1 ? "Replicate $rep" : ""
        scatter!(pl_fit[1], data.tH, data.H;
            label=lbl, color=replicate_colors[rep], alpha=0.7,
            ylabel="Host\n[cell/mL]",
            ylims=(1e2, 1e9), yticks=(ytick_vals1, ytick_labels1),
            markershape=:circle, markersize=8,
            legend=(0.15, 1))
    end
end

# Model for H
for cycle in 1:cycles_sim
    H0_mean = mean(raw_data[rep][cycle].H[1] for rep in replicates)
    V0_mean = mean(raw_data[rep][cycle].V[1] for rep in replicates)
    u0_mean = [H0_mean, V0_mean]

    times = Float64[]
    for rep in replicates
        append!(times, raw_data[rep][cycle].tH)
        append!(times, raw_data[rep][cycle].tV)
    end
    t_c = sort(unique(times))
    t0_c, t1_c = t_c[1], t_c[end]

    sol = solve(
        ODEProblem(model!, u0_mean, (t0_c, t1_c)),
        Rodas5(),
        reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain
    )

    mask_sol = sol.t .<= t_limit
    idx_sol = findall(mask_sol)
    lbl = cycle == 1 ? "Model H" : ""
    if !isempty(idx_sol)
        plot!(pl_fit[1], sol.t[idx_sol], sol[1, idx_sol], label=lbl, color=model_color, lw=4, alpha=0.7)
    end
end

# --- Data & Model for V ---
for rep in replicates
    for cyc_idx in 1:cycles_sim
        data = raw_data[rep][cyc_idx]
        lbl  = cyc_idx == 1 ? "Replicate $rep" : ""
        scatter!(pl_fit[2], data.tV, data.V;
            label=lbl, color=replicate_colors[rep], alpha=0.7,
            ylabel="Virus\n[part/mL]",
            ylims=(1e5, 1e10), yticks=(ytick_vals2, ytick_labels2),
            markershape=:circle, markersize=8,
            legend=(0.15, 0.4))
    end
end

# Model for V
for cycle in 1:cycles_sim
    H0_mean = mean(raw_data[rep][cycle].H[1] for rep in replicates)
    V0_mean = mean(raw_data[rep][cycle].V[1] for rep in replicates)
    u0_mean = [H0_mean, V0_mean]

    times = Float64[]
    for rep in replicates
        append!(times, raw_data[rep][cycle].tH)
        append!(times, raw_data[rep][cycle].tV)
    end
    t_c = sort(unique(times))
    t0_c, t1_c = t_c[1], t_c[end]

    sol = solve(
        ODEProblem(model!, u0_mean, (t0_c, t1_c)),
        Rodas5(),
        reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain
    )

    mask_sol = sol.t .<= t_limit
    idx_sol = findall(mask_sol)
    lbl = cycle == 1 ? "Model V" : ""
    if !isempty(idx_sol)
        plot!(pl_fit[2], sol.t[idx_sol], sol[2, idx_sol], label=lbl, color=model_color, lw=4, alpha=0.7)
    end
end

# --- φ with knots ---
plot!(pl_fit[3], t_plot_limited, ϕ_prime.(t_plot_limited);
      label="ϕ", color=model_color, lw=4, alpha=0.7,
      ylabel="ϕ\n[mL/(part.day)]",
      ylims=(1e-13, 1e-8), yticks=(ytick_vals3, ytick_labels3))

# Add knots to phi subplot (only those up to t_limit)
knots_limited = filter(k -> k <= t_limit, knots)
phi_knot_values = ϕ_prime.(knots_limited)
scatter!(pl_fit[3], knots_limited, phi_knot_values,
    color=:red, markersize=8, markershape=:diamond,
    label="Knots", legend=(0.105, 0.25))

# --- Calcul des bornes temporelles par cycle ---
cycle_tbounds = Vector{Tuple{Float64,Float64}}(undef, cycles_sim)
for cycle in 1:cycles_sim
    t0 = minimum(raw_data[rep][cycle].tH[1] for rep in replicates)
    t1 = maximum(max(raw_data[rep][cycle].tH[end],
                     raw_data[rep][cycle].tV[end]) for rep in replicates)
    cycle_tbounds[cycle] = (t0, t1)
end

# --- Traits pointillés aux dilutions sur tous les subplots ---
dilution_times = [cycle_tbounds[c][2] for c in 1:cycles_sim-1]

for panel in 1:3
    for (i, td) in enumerate(dilution_times)
        vline!(pl_fit[panel], [td];
               color=data_color, lw=2, ls=:dot,
               label = (i == 1 ? "Dilution" : nothing))
    end
end

display(pl_fit)

mkpath(joinpath(@__DIR__, "output"))
fig_path = joinpath(@__DIR__, "output/phi_$(ϕ_combi)_HV_phi.png")
savefig(pl_fit, fig_path)
println("\nFigure saved to $fig_path")
