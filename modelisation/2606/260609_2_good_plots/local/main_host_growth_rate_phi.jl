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


## ===== Colors =====
colors = [
    RGB(0.95, 0.45, 0.45),  # rose
    RGB(0.95, 0.70, 0.30),  # pêche
    RGB(0.40, 0.78, 0.40),  # vert
    RGB(0.35, 0.60, 0.95),  # bleu
    RGB(0.70, 0.45, 0.95),  # lavande
]


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
ytick_vals_mop   = [10.0^i for i in -2:0]
ytick_labels_mop = [L"10^{%$i}" for i in -2:0]
xtick_vals3   = [10.0^i for i in -13:-8]
xtick_labels3 = [L"10^{%$i}" for i in -13:-8]

pl_fit = plot(
    size = (1200, 800),
    left_margin = 15mm,
    right_margin = 10mm,
    top_margin = 10mm,
    bottom_margin = 10mm,
    grid = true,
    xscale = :log10,
    yscale = :log10,
    ylabel = "φ [mL/(part.day)]",
    xlabel = "Host growth rate [1/day]",
    title = "MOP vs φ (φ: $ϕ_combi)",
    legendfontsize = 15,
    guidefontsize = 20,
    tickfontsize = 15,
    titlefontsize = 20,
    legend = :bottomleft
)

# Solve model for each cycle and collect (phi, MOP) pairs
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

    # Calculate φ and MOP at each solution time point
    phi_vals = ϕ_prime.(sol.t)
    growth_rate_vals = r .* (1 .- sol[1,:] ./ K)
    # Filter out null values (zeros, infinities, NaNs)
    valid_idx = findall(x -> x > 0 && isfinite(x), phi_vals) ∩ findall(x -> x > 0 && isfinite(x), growth_rate_vals)
    phi_vals = phi_vals[valid_idx]
    growth_rate_vals = growth_rate_vals[valid_idx]

    plot!(pl_fit, growth_rate_vals, phi_vals, label="Cycle $cycle", color=colors[cycle], lw=3, alpha=0.7,
    xlims=(1e-2, 1e0), xticks=(ytick_vals_mop, ytick_labels_mop),
    ylims=(1e-13, 1e-8), yticks=(xtick_vals3, xtick_labels3))

    scatter!(pl_fit,
         [growth_rate_vals[1]],
         [phi_vals[1]],
         marker=:circle,
         ms=5,
         markerstrokecolor=:transparent,
         markerstrokewidth=0,
         color=colors[cycle],
         label="")

    scatter!(pl_fit,
            [growth_rate_vals[end]],
            [phi_vals[end]],
            marker=:square,
            ms=5,
            markerstrokecolor=:transparent,
            markerstrokewidth=0,
            color=colors[cycle],
            label="")
end

scatter!(pl_fit,
    [0.1],
    [1],
    marker=:circle,
    ms=1,
    markerstrokecolor=:transparent,
    markerstrokewidth=0,
    color=:black,
    label= "Start")

scatter!(pl_fit,
    [5e-13],
    [1],
    marker=:square,
    ms=1,
    markerstrokecolor=:transparent,
    markerstrokewidth=0,
    color=:black,
    label= "End")

display(pl_fit)

mkpath(joinpath(@__DIR__, "output"))
fig_path = joinpath(@__DIR__, "output/host_growth_rate_vs_phi_$(ϕ_combi).png")
savefig(pl_fit, fig_path)
println("\nFigure saved to $fig_path")
