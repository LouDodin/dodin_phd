## =====================================================================
## Burst size — méthode 1 : calcul direct par intervalle de temps
##
## Principe : pour chaque paire de points consécutifs (t_i, t_i+1),
##   burst(t) = ΔV / (-ΔH)        si H diminue (lyse) entre t_i et t_i+1
##
## V est interpolé linéairement sur la grille temporelle de H (tH),
## car les deux fichiers (host/virus) n'ont pas forcément exactement
## les mêmes temps d'échantillonnage.
##
## Suppose que `raw_data` a déjà été construit comme dans votre script
## principal : raw_data[rep][cyc] = (index, tH, H, tV, V)
## =====================================================================

using CSV
using DataFrames
using Plots
using Measures
using Statistics


const replicates  = ["A", "B", "C"]
const cycles_sim  = 5

output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)


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

# Check plot
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
#display(pl_data)




## ----- Interpolation linéaire simple (évite d'ajouter Interpolations.jl) -----
function lin_interp(x_new::Real, x::Vector{Float64}, y::Vector{Float64})
    if x_new <= x[1]
        return y[1]
    elseif x_new >= x[end]
        return y[end]
    end
    i = searchsortedlast(x, x_new)
    i = clamp(i, 1, length(x)-1)
    x0, x1 = x[i], x[i+1]
    y0, y1 = y[i], y[i+1]
    return y0 + (y1 - y0) * (x_new - x0) / (x1 - x0)
end

## ----- Calcul du burst size par intervalle pour un (rep, cycle) donné -----
"""
    compute_burst_direct(tH, H, tV, V; min_relative_decline=0.0)

Retourne un DataFrame avec :
- t_mid       : temps médian de l'intervalle (jour)
- dH          : -ΔH (nombre de cellules "perdues" = lysées sur l'intervalle)
- dV          : ΔV (virus produits net sur l'intervalle)
- burst       : dV/dH, NaN si dH <= 0 (pas de lyse nette mesurable)
- H_start     : H au début de l'intervalle (pour info / pondération)

`min_relative_decline` permet de filtrer les intervalles où la baisse de H
est trop faible par rapport au bruit de mesure (ex : 0.02 = on ignore les
baisses de moins de 2% de H_start, considérées comme du bruit).
"""
function compute_burst_direct(tH::Vector{Float64}, H::Vector{Float64},
                               tV::Vector{Float64}, V::Vector{Float64};
                               min_relative_decline::Float64=0.0)

    n = length(tH)
    V_on_H = [lin_interp(t, tV, V) for t in tH]

    t_mid  = Float64[]
    dH_v   = Float64[]
    dV_v   = Float64[]
    burst  = Float64[]
    H_start = Float64[]

    for i in 1:(n-1)
        ΔH = H[i] - H[i+1]                 # positif si lyse nette
        ΔV = V_on_H[i+1] - V_on_H[i]       # positif si production nette

        push!(t_mid, (tH[i] + tH[i+1]) / 2)
        push!(dH_v, ΔH)
        push!(dV_v, ΔV)
        push!(H_start, H[i])

        rel_decline = ΔH / H[i]
        if ΔH > 0 && rel_decline >= min_relative_decline
            push!(burst, ΔV / ΔH)
        else
            push!(burst, NaN)
        end
    end

    return DataFrame(t_mid=t_mid, dH=dH_v, dV=dV_v, burst=burst, H_start=H_start)
end

## ----- Application à toutes les répliques / cycles -----
## (réutilise `raw_data`, `replicates`, `cycles_sim`, `replicate_colors`,
##  `output_dir` déjà définis dans votre script principal)

burst_results = Dict{String, Vector{DataFrame}}()

for rep in replicates
    dfs = DataFrame[]
    for cyc in 1:length(raw_data[rep])
        data = raw_data[rep][cyc]
        df_b = compute_burst_direct(data.tH, data.H, data.tV, data.V;
                                     min_relative_decline=0.02)
        df_b.cycle = fill(data.index, nrow(df_b))
        push!(dfs, df_b)
    end
    burst_results[rep] = dfs
end


## ----- Visualisation : burst size en fonction du temps -----
pl_burst = plot(size=(700,450), margins=5mm,
                xlabel="Temps (jours)", ylabel="Burst size (virions/cellule lysée)",
                title="Burst size estimé par intervalle (méthode directe)",
                legend=:topright)

for rep in replicates
    for df_b in burst_results[rep]
        valid = .!isnan.(df_b.burst)
        scatter!(pl_burst, df_b.t_mid[valid], df_b.burst[valid];
                 color=replicate_colors[rep],
                 label = df_b.cycle[1]==1 ? "Rep $rep" : "",
                 markersize=4, markerstrokewidth=0)
    end
end

savefig(pl_burst, joinpath(output_dir, "burst_size_direct_vs_time.png"))
display(pl_burst)