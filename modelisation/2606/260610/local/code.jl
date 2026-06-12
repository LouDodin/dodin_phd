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






## ===== Temps commun + moyenne des réplicats =====

function average_over_replicates(replicates, raw_data, t_field::Symbol, val_field::Symbol)
    t_all = sort(unique(vcat([vcat([e[t_field] for e in raw_data[rep]]...) for rep in replicates]...)))
    
    val_all = Vector{Float64}(undef, length(t_all))
    for (i, t) in enumerate(t_all)
        vals = Float64[]
        for rep in replicates
            t_rep   = vcat([e[t_field]   for e in raw_data[rep]]...)
            val_rep = vcat([e[val_field] for e in raw_data[rep]]...)
            idx = findfirst(x -> isapprox(x, t; atol=1e-10), t_rep)
            idx !== nothing && push!(vals, val_rep[idx])
        end
        val_all[i] = mean(vals)
    end
    return t_all, val_all
end

t_H_all, H_all = average_over_replicates(replicates, raw_data, :tH, :H)
t_V_all, V_all = average_over_replicates(replicates, raw_data, :tV, :V)

t_all  = sort(unique(vcat(t_H_all, t_V_all)))  # fixed: vcat instead of two-arg unique
phi_all = ϕ_prime.(t_all)



## ===== Comparaison t_H_all vs t_V_all =====

only_in_H = setdiff(t_H_all, t_V_all)
only_in_V = setdiff(t_V_all, t_H_all)
in_both   = intersect(t_H_all, t_V_all)

println("=== Résumé ===")
println("  t_H_all : $(length(t_H_all)) points")
println("  t_V_all : $(length(t_V_all)) points")
println("  En commun       : $(length(in_both))")
println("  Seulement en H  : $(length(only_in_H))")
println("  Seulement en V  : $(length(only_in_V))")

if !isempty(only_in_H)
    println("\n--- Points présents dans t_H_all mais pas t_V_all ---")
    for t in only_in_H
        println("  t = $t")
    end
end

if !isempty(only_in_V)
    println("\n--- Points présents dans t_V_all mais pas t_H_all ---")
    for t in only_in_V
        println("  t = $t")
    end
end

# Vérification des quasi-doublons (même valeur à la précision flottante près)
println("\n--- Quasi-doublons potentiels (|Δt| < 1e-8) entre H-only et V-only ---")
found = false
for th in only_in_H, tv in only_in_V
    if abs(th - tv) < 1e-8
        println("  t_H=$th  ≈  t_V=$tv  (Δ=$(th-tv))")
        found = true
    end
end
found || println("  Aucun")


## ===== Export CSV =====

CSV.write("modelisation/output/H.csv",
    DataFrame(t = t_H_all, H = H_all))

CSV.write("modelisation/output/V.csv",
    DataFrame(t = t_V_all, V = V_all))

CSV.write("modelisation/output/phi.csv",
    DataFrame(t = t_all, phi = phi_all))

println("Export terminé :")