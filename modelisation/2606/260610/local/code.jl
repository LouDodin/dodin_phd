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






## ===== Moyenne des réplicats par cycle =====

function average_over_replicates_cycle(replicates, raw_data, cyc_idx, t_field::Symbol, val_field::Symbol)
    t_all = sort(unique(vcat([
        raw_data[rep][cyc_idx][t_field]
        for rep in replicates if cyc_idx <= length(raw_data[rep])
    ]...)))

    val_all = Vector{Float64}(undef, length(t_all))
    for (i, t) in enumerate(t_all)
        vals = Float64[]
        for rep in replicates
            cyc_idx > length(raw_data[rep]) && continue
            t_rep   = raw_data[rep][cyc_idx][t_field]
            val_rep = raw_data[rep][cyc_idx][val_field]
            idx = findfirst(x -> isapprox(x, t; atol=1e-10), t_rep)
            idx !== nothing && push!(vals, val_rep[idx])
        end
        val_all[i] = mean(vals)
    end
    return t_all, val_all
end


## ===== Export CSV par cycle =====

mkpath("modelisation/output")

for cyc_idx in 1:cycles_sim
    # Vérifie qu'au moins un réplicat a ce cycle
    any(cyc_idx <= length(raw_data[rep]) for rep in replicates) || continue

    t_H, H_cyc = average_over_replicates_cycle(replicates, raw_data, cyc_idx, :tH, :H)
    t_V, V_cyc = average_over_replicates_cycle(replicates, raw_data, cyc_idx, :tV, :V)

    t_phi = sort(unique(vcat(t_H, t_V)))
    phi_cyc = ϕ_prime.(t_phi)

    CSV.write("modelisation/output/H_cycle$(cyc_idx).csv",
        DataFrame(t = t_H, H = H_cyc))

    CSV.write("modelisation/output/V_cycle$(cyc_idx).csv",
        DataFrame(t = t_V, V = V_cyc))

    CSV.write("modelisation/output/phi_cycle$(cyc_idx).csv",
        DataFrame(t = t_phi, phi = phi_cyc))

    println("Cycle $cyc_idx exporté.")
end

println("Export terminé.")