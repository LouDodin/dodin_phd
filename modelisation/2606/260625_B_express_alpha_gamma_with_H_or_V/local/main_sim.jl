## ===== Packages =====
using DifferentialEquations
using OrdinaryDiffEqRosenbrock
using CSV
using DataFrames
using Statistics
using Plots
using Measures
using LaTeXStrings
using Colors


## ===== Choices =====
const MODEL_FILE = joinpath(@__DIR__, "models/SR_RS_H.jl")
const replicates  = ["A", "B", "C"]
const cycles_sim  = 5


## ===== Utils =====
include(MODEL_FILE)
using .ModelDef

output_dir = joinpath(@__DIR__, "output/$(ModelDef.MODEL_NAME)")
mkpath(output_dir)

const isoutofdomain = (u, p, t) -> any(x -> x < 0, u)

replicate_colors = Dict(
    "A" => RGB(0.6, 0.8, 1.0),
    "B" => RGB(31/255, 119/255, 180/255),
    "C" => RGB(0.0, 0.3, 0.7)
)

model_color = RGB(255/255, 127/255, 14/255)
cycle_color = RGB(31/255, 119/255, 180/255)


## ===== Parameters (defined manually) =====
p_best = [
    0.2775277333185582, # eps_H/V / c
    0.00778760350728357, # K_H/V / p
    7.478216159796185e-9 # phi
]


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


## ===== Phi_ref =====
const poly_file = "modelisation/2606/260604/genotoul/output/nint_3_2_2_3_2/polynomial.txt"

lines = readlines(poly_file)

function extract_float(s::AbstractString)::Float64
    m = match(r"([+-]?\s*[0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)", strip(s))
    m === nothing && error("Cannot extract float from: \"$s\"")
    parse(Float64, replace(m.captures[1], " " => ""))
end

ϕ_intervals = Vector{NTuple{6,Float64}}()
global i = 1
while i <= length(lines)
    m_iv = match(r"Interval\s+\[([0-9eE+\-.]+),\s*([0-9eE+\-.]+)\]\s+days:", strip(lines[i]))
    if m_iv !== nothing
        if i + 4 > length(lines)
            global i += 1; continue
        end
        t0 = parse(Float64, m_iv.captures[1])
        t1 = parse(Float64, m_iv.captures[2])
        a  = extract_float(split(lines[i+1], "=")[end])
        b  = extract_float(split(lines[i+2], "*")[1])
        c  = extract_float(split(lines[i+3], "*")[1])
        d  = extract_float(split(lines[i+4], "*")[1])
        push!(ϕ_intervals, (t0, t1, a, b, c, d))
        i += 5; continue
    end
    global i += 1
end

const tmin  = ϕ_intervals[1][1]
const tmax  = ϕ_intervals[end][2]

function ϕ_ref(t)
    tc  = clamp(t, tmin, tmax)
    idx = length(ϕ_intervals)
    for k in eachindex(ϕ_intervals)
        if tc <= ϕ_intervals[k][2]; idx = k; break; end
    end
    t0, _, a, b, c, d = ϕ_intervals[idx]
    dt = tc - t0
    exp(a + b*dt + c*dt^2 + d*dt^3)
end

const N_GRID = 10_000
const t_grid = collect(range(tmin, tmax; length=N_GRID))
const ϕ_grid = ϕ_ref.(t_grid)

knots = Float64[]
for interval in ϕ_intervals
    push!(knots, interval[1])
end
push!(knots, ϕ_intervals[end][2])
unique!(sort!(knots))


## ===== Simulation & Plot =====
ytick_vals1   = [10.0^i for i in 2:9]
ytick_labels1 = [L"10^{%$i}" for i in 2:9]
ytick_vals2   = [10.0^i for i in 5:10]
ytick_labels2 = [L"10^{%$i}" for i in 5:10]
ytick_vals3   = [10.0^i for i in -14:-7]
ytick_labels3 = [L"10^{%$i}" for i in -14:-7]
ytick_vals4   = [10.0^i for i in -8:-3]
ytick_labels4 = [L"10^{%$i}" for i in -8:-3]

param_str = join(
    ["$(ModelDef.FITTED_PARAMS[i].name)=$(round(p_best[i], sigdigits=2))"
     for i in eachindex(ModelDef.FITTED_PARAMS)],
    "  "
)

pl_fit = plot(
    layout = (3, 2),
    size = (2500, 2000),
    left_margin = 20mm,
    right_margin = 10mm,
    top_margin = 10mm,
    bottom_margin = 15mm,
    grid = true,
    yscale = :log10,
    xlims = (0, 67),
    ytickfontsize = 25,
    legendfontsize = 18,
    guidefontsize = 23,
    xtickfontsize = 23,
    titlefontsize = 23,
    xlabel = "Time [days]",
    legend = :bottomright,
    plot_title = "$(ModelDef.MODEL_NAME)\n$(param_str)",
    plot_titlefontsize = 28
)

# --- Data ---
for rep in replicates
    for cyc_idx in 1:cycles_sim
        data = raw_data[rep][cyc_idx]
        lbl  = cyc_idx == 1 ? "Replicate $rep" : ""
        scatter!(pl_fit[1], data.tH, data.H;
            label=lbl, color=replicate_colors[rep], alpha=0.7,
            ylabel="Host\n[cell/mL]",
            ylims=(1e2, 1e9), yticks=(ytick_vals1, ytick_labels1),
            markershape=:circle, markersize=13,
            legend=(0.15, 1))
        scatter!(pl_fit[3], data.tV, data.V;
            label=lbl, color=replicate_colors[rep], alpha=0.7,
            ylabel="Virus\n[part/mL]",
            ylims=(1e5, 1e10), yticks=(ytick_vals2, ytick_labels2),
            markershape=:circle, markersize=13,
            legend=(0.15, 0.3))
    end
end

plot!(pl_fit[5], t_grid, ϕ_grid;
      label="ϕ_ref", color=replicate_colors["B"], lw=8, alpha=0.7,
      ylabel="ϕ\n[mL/(part.day)]",
      ylims=(1e-14, 1e-7), yticks=(ytick_vals3, ytick_labels3))

phi_knot_values = ϕ_ref.(knots)
scatter!(pl_fit[5], knots, phi_knot_values, color=replicate_colors["C"],
      markersize=13, markershape=:diamond,
      label="Knots", legend=:bottomright)

# --- Simulation ---
global prop_S_plot = 1.0

for cycle in 1:cycles_sim
    H0 = mean(raw_data[rep][cycle].H[1] for rep in replicates)
    V0 = mean(raw_data[rep][cycle].V[1] for rep in replicates)
    u0 = ModelDef.INITIAL_CONDITION(H0, V0, prop_S_plot)

    t0 = minimum(raw_data[rep][cycle].tH[1]   for rep in replicates)
    t1 = maximum(max(raw_data[rep][cycle].tH[end], raw_data[rep][cycle].tV[end]) for rep in replicates)

    sol = solve(ODEProblem(ModelDef.ODE_MODEL!, u0, (t0, t1), p_best),
                Rodas5(), reltol=1e-6, abstol=1e-6,
                isoutofdomain=isoutofdomain)

    u_end = sol.u[end]
    H_end = u_end[1] + u_end[2]
    global prop_S_plot = H_end > 0 ? clamp(u_end[1] / H_end, 0.0, 1.0) : 1.0

    S = [u[1] for u in sol.u]
    R = [u[2] for u in sol.u]
    V = [u[3] for u in sol.u]
    H = S .+ R
    ϕ_equiv = [ModelDef.PHI_EQUIV(sol(t), p_best) for t in sol.t]
    if ModelDef.MODEL_NAME == "SR_RS_H"
        alpha = p_best[1] .* H ./ (H .+ p_best[2])
        gamma = p_best[1] .* p_best[2] ./ (H .+ p_best[2])
    elseif ModelDef.MODEL_NAME == "SR_RS_V"
        alpha = p_best[1] .* V ./ (V .+ p_best[2])
        gamma = p_best[1] .* p_best[2] ./ (V .+ p_best[2])
    elseif ModelDef.MODEL_NAME == "SR_RS"
        alpha = fill(p_best[1] * p_best[2], length(H))
        gamma = fill(p_best[1] * (1-p_best[2]), length(H))
    end

    plot!(pl_fit[1], sol.t, H, label=cycle==1 ? "H model" : "", color=model_color, lw=8, alpha=0.7)
    plot!(pl_fit[1], sol.t, S, label=cycle==1 ? "S model" : "", color=:red, lw=6, ls=:dash, alpha=0.7)
    plot!(pl_fit[1], sol.t, R, label=cycle==1 ? "R model" : "", color=:green, lw=6, ls=:dash, alpha=0.7)
    plot!(pl_fit[2], sol.t, alpha, label=cycle==1 ? "Alpha model" : "", color=model_color, lw=8, alpha=0.7, ylabel="Alpha\n[1/day]", ylims=(1e-8, 1e-3), yticks=(ytick_vals4, ytick_labels4))
    plot!(pl_fit[3], sol.t, V,     label=cycle==1 ? "V model" : "",     color=model_color, lw=8, alpha=0.7)
    plot!(pl_fit[4], sol.t, gamma, label=cycle==1 ? "Gamma model" : "", color=model_color, lw=8, alpha=0.7, ylabel="Gamma\n[1/day]", yscale=:identity)
    plot!(pl_fit[5], sol.t, ϕ_equiv, label=cycle==1 ? "ϕ_equiv model" : "", color=model_color, lw=8, alpha=0.7)
end

# --- Dilutions ---
cycle_tbounds = Vector{Tuple{Float64,Float64}}(undef, cycles_sim)
for cycle in 1:cycles_sim
    t0 = minimum(raw_data[rep][cycle].tH[1] for rep in replicates)
    t1 = maximum(max(raw_data[rep][cycle].tH[end],
                     raw_data[rep][cycle].tV[end]) for rep in replicates)
    cycle_tbounds[cycle] = (t0, t1)
end

dilution_times = [cycle_tbounds[c][2] for c in 1:cycles_sim-1]

for panel in 1:5
    for td in dilution_times
        vline!(pl_fit[panel], [td]; color=cycle_color, lw=4, ls=:dot, label=nothing)
    end
end

display(pl_fit)

fig_path = "$(output_dir)/$(ModelDef.MODEL_NAME)_sim.png"
savefig(pl_fit, fig_path)
println("\nFigure saved to $fig_path")