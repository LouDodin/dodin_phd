## ===== Packages =====
using CSV
using DataFrames
using Statistics
using Plots
using Measures
using DifferentialEquations
using Interpolations
using DataInterpolations

## ===== Choices =====
ϕ_fix = 1e-6
const ϕ_combi     = "3_2_2_3_2"
const replicates  = ["A", "B", "C"]
const cycles_sim  = 5


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






## ===== Temps commun + moyenne des réplicats =====

t_H_all = sort(unique(vcat([vcat([e.tH for e in raw_data[rep]]...) for rep in replicates]...)))

println("Nombre de points de temps communs : $(length(t_H_all))")

H_all = Vector{Float64}(undef, length(t_H_all))

for (i, t) in enumerate(t_H_all)
    vals = Float64[]
    for rep in replicates
        t_all = vcat([e.tH for e in raw_data[rep]]...)
        H_all_rep = vcat([e.H  for e in raw_data[rep]]...)
        idx = findfirst(t_all .== t)
        idx !== nothing && push!(vals, H_all_rep[idx])
    end
    H_all[i] = mean(vals)
end

println("Longueur H_all : $(length(H_all))")






## --- φ'(t) polynomial ---
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

# Check φ'(t)
tmin   = ϕ_intervals[1][1]
tmax   = ϕ_intervals[end][2]
t_plot = range(tmin, tmax, length=5000)
display(plot(t_plot, ϕ_prime.(t_plot), lw=2, xlabel="time (days)", ylabel="ϕ'(t)", title="ϕ'(t) polynomial", yaxis=:log))








## ===== Calculs =====
## ===== Calculs =====
ϕ_prime_all = ϕ_prime.(t_H_all)

R_all = H_all .* (1 .- ϕ_prime_all ./ ϕ_fix)

S_all = H_all .- R_all

r = 0.574619342477644
K = 6.675449070379925e7
β = 144.0
δ = 0.02

function itp_S(t)
    tc = clamp(t, t_H_all[1], t_H_all[end])
    i = clamp(searchsortedlast(t_H_all, tc), 1, length(t_H_all) - 1)
    t0, t1 = t_H_all[i], t_H_all[i+1]
    s0, s1 = S_all[i], S_all[i+1]
    s0 + (s1 - s0) * (tc - t0) / (t1 - t0)
end

p = (β, δ, ϕ_fix)





## ===== Vi =====
cycle_bounds = [(raw_data["A"][cyc].tH[1], raw_data["A"][cyc].tH[end]) for cyc in 1:cycles_sim]

function ode_Vi!(dY, Y, p, t)
    β, δ, ϕ = p
    dY[1] = β * ϕ * itp_S(t) * Y[1] - δ * Y[1]
end

Vi_sim = Float64[]
t_sim  = Float64[]

for cyc in 1:cycles_sim
    t_start, t_end = cycle_bounds[cyc]

    mask   = t_start .<= t_H_all .<= t_end
    t_cyc  = t_H_all[mask]
    isempty(t_cyc) && continue

    Vi0_cyc = mean(
        raw_data[rep][cyc].V[1]
        for rep in replicates
        if cyc <= length(raw_data[rep])
    )

    prob2 = ODEProblem(ode_Vi!, [Vi0_cyc], (t_cyc[1], t_cyc[end]), p)
    sol2  = solve(prob2, Tsit5(), saveat=t_cyc)

    append!(t_sim,  t_cyc)
    append!(Vi_sim, sol2[1, :])
end





## ===== Plot =====

pl_HVR = plot(
    layout = (5, 1),
    size = (1800, 2200),
    left_margin = 15mm,
    right_margin = 10mm,
    top_margin = 8mm,
    bottom_margin = 10mm,
    grid = true,
    legend = :bottomright,
    xlims = (0, 67),
    xlabel = "Time (days)",
    xtickfontsize = 18,
    ytickfontsize = 18,
    guidefontsize = 20,
    legendfontsize = 14,
    titlefontsize = 20,
)

# =========================
# 1) HOST H
# =========================
for rep in replicates
    for cyc in 1:cycles_sim
        data = raw_data[rep][cyc]
        scatter!(pl_HVR[1], data.tH, data.H;
            color = replicate_colors[rep],
            alpha = 0.7,
            markersize = 6,
            yscale = :log10,
            ylabel = "Host (cell/mL)",
            ylims = (1e2, 1e8),
            label = cyc == 1 ? "Rep $rep" : "",
        )
    end
end

scatter!(pl_HVR[1], t_H_all, H_all;
        color = :red,
        alpha = 0.7,
        markersize = 6,
        label = "H_all"
    )


# =========================
# 2) VIRUS V
# =========================
for rep in replicates
    for cyc in 1:cycles_sim
        data = raw_data[rep][cyc]
        scatter!(pl_HVR[2], data.tV, data.V;
            color = replicate_colors[rep],
            alpha = 0.7,
            markersize = 6,
            yscale = :log10,
            ylabel = "Virus (part/mL)",
            ylims = (1e3, 1e10),
            label = cyc == 1 ? "Rep $rep" : "",
        )
    end
end

plot!(pl_HVR[2], t_sim, Vi_sim;
    lw    = 2,
    color = :black,
    label = "Vi simulé",
    yscale = :log10,
)

# =========================
# 3) PHI
# =========================
plot!(pl_HVR[3], t_plot, ϕ_prime.(t_plot);
    lw = 2,
    color = :red,
    label = "φ'",
    ylabel = "φ'",
    yscale = :log10,
)

# =========================
# 4) R = H (1 - φ'/φ_fix)
# =========================

scatter!(pl_HVR[4], t_H_all, R_all;
        color = :red,
        alpha = 0.7,
        markersize = 6,
        ylabel = "R",
        label = "R",
        yscale = all(>(0), R_all) ? :log10 : :identity
    )


# =========================
# 5) S = itp_S(t)
# =========================
t_itp_plot = range(t_H_all[1], t_H_all[end], length=5000)

"""
plot!(pl_HVR[5], t_itp_plot, itp_S.(t_itp_plot);
    lw     = 2,
    color  = :green,
    label  = "itp_S",
)
"""

scatter!(pl_HVR[5], t_H_all, S_all;
    color      = :darkgreen,
    alpha      = 0.7,
    markersize = 6,
    label      = "S_all",
    ylabel = "S (cell/mL)",
    xlabel = "Time (days)",
    yscale = :log10,
)
# =========================
# display + save
# =========================
display(pl_HVR)

mkpath(joinpath(@__DIR__, "output"))
savefig(pl_HVR, joinpath(@__DIR__, "output/$(ϕ_fix).png"))