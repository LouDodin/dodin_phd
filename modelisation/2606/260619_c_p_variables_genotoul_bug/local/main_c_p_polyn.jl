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
const cp_combi    = "1_1_1_1_1"   # combinaison de c et p
const phi_ref_combi = "3_2_2_3_2" # combinaison du φ_ref
const replicates  = ["A", "B", "C"]
const cycles_fit  = 5
const cycles_sim  = 5
const n_runs      = 1

# Limit display to a specific knot (nothing = no limit, or specify knot index)
const knot_limit = nothing  # Set to knot index (e.g., 3) to limit display


## ===== Colors =====
color_A = RGB(0.6, 0.8, 1.0)
color_B = RGB(31/255, 119/255, 180/255)
color_C = RGB(0.0, 0.3, 0.7)
model_color     = RGB(255/255, 127/255, 14/255)
color_phi_ref   = RGB(31/255, 119/255, 180/255)
data_color      = RGB(31/255, 119/255, 180/255)

replicate_colors = Dict(
    "A" => color_A,
    "B" => color_B,
    "C" => color_C
)


## ===== Inputs =====

# --- H & V data ---
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


# --- Polynôme c(t) et p(t) + ϕ constant ---
const poly_file = joinpath(@__DIR__,
    "../../260619_c_p_variables/genotoul/output/nint_$(cp_combi)/polynomial.txt")

"""
Parse un fichier polynomial.txt produit par fit_c_p.jl.
Retourne (c_func, p_func, ϕ_val, metadata, c_intervals, p_intervals).
Chaque intervalle est un NTuple{6,Float64} : (t0, t1, a, b, c_coef, d).
La fonction évalue exp(a + b*dt + c*dt² + d*dt³) avec dt = t - t0.
"""
function parse_cp_polynome(filepath::String)
    isfile(filepath) || error("Polynomial file not found: $filepath")
    lines = readlines(filepath)

    # ── METADATA ──
    metadata = Dict{String,Any}()
    ϕ_val = NaN
    for line in lines
        # ϕ en valeur naturelle : "ϕ      = 1.234e-10  (mL/cell/day)"
        m = match(r"^ϕ\s*=\s*([\d.e+\-]+)", strip(line))
        if m !== nothing
            ϕ_val = parse(Float64, m.captures[1])
        end
        m2 = match(r"^(total_knots|n_runs|best_fitness|best_seed)\s*=\s*(.+)$", strip(line))
        m2 === nothing && continue
        key = m2.captures[1]
        val = strip(m2.captures[2])
        metadata[key] = key in ("total_knots","n_runs","best_seed") ?
                         parse(Int, val) : parse(Float64, val)
    end
    isnan(ϕ_val) && @warn "ϕ not found in $filepath — defaulting to NaN"

    # ── Extracteur de flottant robuste ──
    function extract_float(s::AbstractString)::Float64
        m = match(r"([+-]?\s*[0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)", strip(s))
        m === nothing && error("Cannot extract float from: \"$s\"")
        parse(Float64, replace(m.captures[1], " " => ""))
    end

    # ── Parseur générique d'une section PIECEWISE POLYNOMIAL ──
    # Cherche les blocs "Interval [t0, t1] days:" appartenant à la section `label`.
    function parse_section(label::String)::Vector{NTuple{6,Float64}}
        intervals = NTuple{6,Float64}[]
        in_section = false
        i = 1
        while i <= length(lines)
            l = strip(lines[i])
            # Détection de la section
            if occursin("PIECEWISE POLYNOMIAL OF log($label(t))", l)
                in_section = true
                i += 1; continue
            end
            # Fin de section : une autre section ========== commence
            if in_section && occursin("==========", l) && !isempty(l)
                break
            end
            if in_section
                m_iv = match(r"Interval\s+\[([0-9eE+\-.]+),\s*([0-9eE+\-.]+)\]\s+days:", l)
                if m_iv !== nothing
                    if i + 4 > length(lines)
                        @warn "Incomplete interval block at line $i — skipping"
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
            end
            i += 1
        end
        isempty(intervals) && error("No intervals found for $label in $filepath")
        return intervals
    end

    c_intervals = parse_section("c")
    p_intervals = parse_section("p")

    # ── Constructeur de fonction évaluatrice ──
    function make_func(intervals::Vector{NTuple{6,Float64}})
        t_lo = intervals[1][1]
        t_hi = intervals[end][2]
        return function(t)
            tc  = clamp(t, t_lo, t_hi)
            idx = length(intervals)
            for k in eachindex(intervals)
                if tc <= intervals[k][2]; idx = k; break; end
            end
            t0, _, a, b, c_coef, d = intervals[idx]
            dt = tc - t0
            exp(a + b*dt + c_coef*dt^2 + d*dt^3)
        end
    end

    c_func = make_func(c_intervals)
    p_func = make_func(p_intervals)

    return c_func, p_func, ϕ_val, metadata, c_intervals, p_intervals
end

const c_func, p_func, ϕ_val, cp_meta, c_intervals, p_intervals =
    parse_cp_polynome(poly_file)

println("ϕ (constant) = $ϕ_val")
println("Metadata: $cp_meta")


# --- ϕ_ref : polynôme de référence (ancien modèle phi variable) ---
const phi_ref_file = joinpath(@__DIR__,
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
            if i + 4 > length(lines)
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

# ── Grille temporelle commune (nœuds de c et p sont identiques par construction) ──
tmin = c_intervals[1][1]
tmax = c_intervals[end][2]

# Interpolations linéaires sur grille fine pour l'ODE (thread-safe, rapide)
const N_GRID   = 10_000
const t_grid   = collect(range(tmin, tmax; length=N_GRID))
const c_grid   = c_func.(t_grid)
const p_grid   = p_func.(t_grid)
const c_interp = LinearInterpolation(c_grid, t_grid)
const p_interp = LinearInterpolation(p_grid, t_grid)

# ── Nœuds (communs à c et p) ──
knots = Float64[]
for iv in c_intervals; push!(knots, iv[1]); end
push!(knots, c_intervals[end][2])
unique!(sort!(knots))

# ── Limite d'affichage selon knot_limit ──
t_limit = if knot_limit === nothing
    tmax
else
    (knot_limit < 1 || knot_limit > length(knots)) &&
        error("knot_limit must be between 1 and $(length(knots))")
    knots[knot_limit]
end
t_plot         = range(tmin, tmax,           length=5000)
t_plot_limited = range(tmin, min(t_limit, tmax), length=5000)


## ===== Constants =====
const r = 0.574619342477644
const K = 6.675449070379925e7
const β = 144.0
const δ = 0.02


## ===== Model SR_RS =====
# États : Y = [S, R, Vi]
# c(t) : taux de changement de sensibilité/résistance
# p(t) : proportion de S qui devient R
# ϕ    : constante d'adsorption

function model!(dY, Y, params, t)
    c_itp, p_itp, ϕ = params
    S, R, Vi = Y[1], Y[2], Y[3]
    H = S + R
    ct = c_itp(clamp(t, t_grid[1], t_grid[end]))
    pt = p_itp(clamp(t, t_grid[1], t_grid[end]))
    dY[1] = r*S*(1 - H/K) - ϕ*S*Vi - ct*pt*S + ct*(1 - pt)*R
    dY[2] = r*R*(1 - H/K) + ct*pt*S - ct*(1 - pt)*R
    dY[3] = β*ϕ*S*Vi - δ*Vi
end

isoutofdomain(u, p, t) = any(x -> x < 0, u)

ode_params = (c_interp, p_interp, ϕ_val)


## ===== Simulation & Plot =====

ytick_vals1   = [10.0^i for i in 2:9]
ytick_labels1 = [L"10^{%$i}" for i in 2:9]
ytick_vals2   = [10.0^i for i in 5:10]
ytick_labels2 = [L"10^{%$i}" for i in 5:10]
ytick_vals3   = [10.0^i for i in -17:-1]
ytick_labels3 = [L"10^{%$i}" for i in -17:-1]
ytick_vals4   = [10.0^i for i in -12:1]
ytick_labels4 = [L"10^{%$i}" for i in -12:1]
ytick_vals5   = [10.0^i for i in -14:-7]
ytick_labels5 = [L"10^{%$i}" for i in -14:-7]

title_suffix = knot_limit === nothing ? "" : " [Knot limit: $knot_limit]"

pl_fit = plot(
    layout = (5, 1),
    size = (1800, 2500),
    left_margin  = 15mm,
    right_margin = 10mm,
    top_margin   = 5mm,
    bottom_margin = 10mm,
    grid = true,
    yscale = :log10,
    xlims = (0, 67),
    ytickfontsize  = 22,
    legendfontsize = 15,
    guidefontsize  = 20,
    xtickfontsize  = 20,
    titlefontsize  = 20,
    xlabel = "Time [days]",
    legend = :bottomright,
    plot_title = "Model SR_RS with variable c(t), p(t) (nint: $cp_combi)" * title_suffix,
    plot_titlefontsize = 25
)

# ── Données scatter (subplots 1 & 2) ────────────────────────────────────────
for rep in replicates
    for cyc_idx in 1:cycles_sim
        data = raw_data[rep][cyc_idx]
        lbl  = cyc_idx == 1 ? "Replicate $rep" : ""
        scatter!(pl_fit[1], data.tH, data.H;
            label=lbl, color=replicate_colors[rep], alpha=0.7,
            ylabel="Host H=S+R\n[cell/mL]",
            ylims=(1e2, 1e9), yticks=(ytick_vals1, ytick_labels1),
            markershape=:circle, markersize=8,
            legend=(0.15, 1))
        scatter!(pl_fit[2], data.tV, data.V;
            label=lbl, color=replicate_colors[rep], alpha=0.7,
            ylabel="Virus Vi\n[part/mL]",
            ylims=(1e5, 1e10), yticks=(ytick_vals2, ytick_labels2),
            markershape=:circle, markersize=8,
            legend=(0.15, 0.4))
    end
end

# ── Résolution ODE par cycle (une seule passe, solutions mémorisées) ──────────
cycle_solutions = Vector{Any}(undef, cycles_sim)

for cycle in 1:cycles_sim
    H0_mean = mean(raw_data[rep][cycle].H[1] for rep in replicates)
    V0_mean = mean(raw_data[rep][cycle].V[1] for rep in replicates)
    u0 = [H0_mean, 0.0, V0_mean]   # S=H0, R=0, Vi=V0

    times = Float64[]
    for rep in replicates
        append!(times, raw_data[rep][cycle].tH)
        append!(times, raw_data[rep][cycle].tV)
    end
    t_c = sort(unique(times))
    t0_c, t1_c = t_c[1], t_c[end]

    sol = solve(
        ODEProblem(model!, u0, (t0_c, t1_c), ode_params),
        Rodas5(),
        reltol=1e-6, abstol=1e-6,
        isoutofdomain=isoutofdomain
    )
    cycle_solutions[cycle] = sol
end

# ── Subplot 1 : H = S + R ────────────────────────────────────────────────────
for cycle in 1:cycles_sim
    sol  = cycle_solutions[cycle]
    mask = sol.t .<= t_limit
    idx  = findall(mask)
    lbl  = cycle == 1 ? "Model H" : ""
    if !isempty(idx)
        plot!(pl_fit[1], sol.t[idx], sol[1, idx] .+ sol[2, idx];
              label=lbl, color=model_color, lw=4, alpha=0.7)
    end
end

# ── Subplot 2 : Virus Vi ─────────────────────────────────────────────────────
for cycle in 1:cycles_sim
    sol  = cycle_solutions[cycle]
    mask = sol.t .<= t_limit
    idx  = findall(mask)
    lbl  = cycle == 1 ? "Model Vi" : ""
    if !isempty(idx)
        plot!(pl_fit[2], sol.t[idx], sol[3, idx];
              label=lbl, color=model_color, lw=4, alpha=0.7)
    end
end

# ── Subplot 3 : α(t) ────────────────────────────────────────────────────────
knots_limited = filter(k -> k <= t_limit, knots)

plot!(pl_fit[3], t_plot_limited, c_func.(t_plot_limited).*p_func.(t_plot_limited);
      label="α(t)", color=model_color, lw=4, alpha=0.8,
      ylabel="α(t)\n[1/day]",
      ylims=(1e-17, 1e-1),
      yticks=(ytick_vals3, ytick_labels3))
scatter!(pl_fit[3], knots_limited, c_func.(knots_limited).*p_func.(knots_limited);
    color=:red, markersize=8, markershape=:diamond,
    label="Knots", legend=(0.105, 0.25))

# ── Subplot 4 : p(t) ────────────────────────────────────────────────────────
plot!(pl_fit[4], t_plot_limited, c_func.(t_plot_limited).*(1 .- p_func.(t_plot_limited));
      label="γ(t)", color=model_color, lw=4, alpha=0.8,
      ylabel="γ(t)\n[dimensionless]",
      ylims=(1e-12, 1e1),
      yticks=(ytick_vals4, ytick_labels4))
scatter!(pl_fit[4], knots_limited, c_func.(knots_limited).*(1 .- p_func.(knots_limited));
    color=:red, markersize=8, markershape=:diamond,
    label="Knots", legend=(0.105, 0.25))

# ── Subplot 5 : ϕ_ref vs ϕ_equiv ────────────────────────────────────────────
# ϕ_ref  : polynôme de référence (ancien modèle phi variable, depuis phi_ref_file)
# ϕ_equiv = ϕ * S / H  : phi effectif du modèle SR_RS
#           (ϕ constant × fraction de sensibles → phi apparent vu du virus)

t_phi_ref_limited = range(tmin, min(t_limit, tmax), length=5000)
plot!(pl_fit[5], t_phi_ref_limited, ϕ_ref_func.(t_phi_ref_limited);
      label="ϕ_ref", color=color_phi_ref, lw=6, alpha=0.7,
      ylabel="ϕ\n[mL/(part.day)]",
      ylims=(1e-14, 1e-7), yticks=(ytick_vals5, ytick_labels5))

knots_ref_limited = filter(k -> k <= t_limit, knots_ref)
scatter!(pl_fit[5], knots_ref_limited, ϕ_ref_func.(knots_ref_limited);
    color=color_phi_ref, markersize=10, markershape=:diamond,
    label="Knots ϕ_ref", legend=:bottomright)

for cycle in 1:cycles_sim
    sol  = cycle_solutions[cycle]
    mask = sol.t .<= t_limit
    idx  = findall(mask)
    lbl  = cycle == 1 ? "ϕ_equiv = ϕ·S/H" : ""
    if !isempty(idx)
        S_vec     = sol[1, idx]
        R_vec     = sol[2, idx]
        H_vec     = S_vec .+ R_vec
        # ϕ_equiv = ϕ * S/H  (→ 0 quand R domine, → ϕ quand S domine)
        phi_equiv = ϕ_val .* S_vec ./ max.(H_vec, 1e-30)
        plot!(pl_fit[5], sol.t[idx], phi_equiv;
              label=lbl, color=model_color, lw=4, alpha=0.7)
    end
end

# ── Traits pointillés aux dilutions ─────────────────────────────────────────
cycle_tbounds = Vector{Tuple{Float64,Float64}}(undef, cycles_sim)
for cycle in 1:cycles_sim
    t0 = minimum(raw_data[rep][cycle].tH[1] for rep in replicates)
    t1 = maximum(max(raw_data[rep][cycle].tH[end],
                     raw_data[rep][cycle].tV[end]) for rep in replicates)
    cycle_tbounds[cycle] = (t0, t1)
end

dilution_times = [cycle_tbounds[c][2] for c in 1:cycles_sim-1]

for panel in 1:5
    for (i, td) in enumerate(dilution_times)
        vline!(pl_fit[panel], [td];
               color=data_color, lw=2, ls=:dot,
               label = (i == 1 ? "Dilution" : nothing))
    end
end

display(pl_fit)

mkpath(joinpath(@__DIR__, "output"))
fig_path = joinpath(@__DIR__, "output/$(cp_combi).png")
savefig(pl_fit, fig_path)
println("\nFigure saved to $fig_path")