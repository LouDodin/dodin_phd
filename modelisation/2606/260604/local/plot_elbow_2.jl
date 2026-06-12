using CSV
using DataFrames
using Plots
using StatsPlots
using Statistics
using Measures
using Printf
using Colors

# ─────────────────────────────────────────────────────────────────────────────
# plot_elbow.jl
#
# Pour chaque cycle, trace :
#   - strip plot de toutes les combinaisons (points avec jitter)
#   - courbe de la médiane et du minimum par valeur de n_interior
#   - effet marginal moyen (Δ%) contrôlé sur les autres cycles
#
# Logique : 1 seed par combinaison → pas de boxplot, on appariera les
# combinaisons qui diffèrent uniquement par n_interior du cycle courant.
#
# Usage :
#   julia plot_elbow.jl
#   julia plot_elbow.jl /chemin/summary_fits.csv
# ─────────────────────────────────────────────────────────────────────────────

# ── Configuration ─────────────────────────────────────────────────────────────

const CSV_PATH  = length(ARGS) > 0 ? ARGS[1] :
                  joinpath(@__DIR__, "../genotoul/summary_fits_12345_seed1.csv")
const N_CYCLES  = 5
const OUT_PATH  = joinpath(@__DIR__, "output/elbow_2.png")

# ── Chargement ────────────────────────────────────────────────────────────────

println("Loading: $CSV_PATH")
df = CSV.read(CSV_PATH, DataFrame)
println("$(nrow(df)) rows loaded")

# ── Parsing de la colonne combination → n_int_c1 … n_int_c5 ──────────────────

function parse_combination(s::AbstractString)::Vector{Int}
    parse.(Int, split(s, "_"))
end

combs = parse_combination.(df.combination)
for cyc in 1:N_CYCLES
    df[!, Symbol("n_int_c$(cyc)")] = getindex.(combs, cyc)
end

# ── Helpers ───────────────────────────────────────────────────────────────────

"""Colonnes des autres cycles (pour l'appariement marginal)."""
other_cols(cyc::Int) = [Symbol("n_int_c$c") for c in 1:N_CYCLES if c != cyc]

"""
    marginal_delta_pct(df, cyc, n)

Calcule l'effet marginal moyen (en %) de passer de n_interior=n à n+1
sur le cycle `cyc`, en appariant les combinaisons identiques sur tous
les autres cycles.

Retourne (Δ_abs, Δ_pct, n_pairs) ou `nothing` si aucune paire trouvée.
"""
function marginal_delta_pct(df::DataFrame, cyc::Int, n::Int)
    col   = Symbol("n_int_c$(cyc)")
    ocols = other_cols(cyc)

    sub_n   = df[df[!, col] .== n,   vcat(ocols, [:fitness])]
    sub_np1 = df[df[!, col] .== n+1, vcat(ocols, [:fitness])]

    isempty(sub_n) || isempty(sub_np1) && return nothing

    joined = innerjoin(sub_n, sub_np1, on = ocols,
                       makeunique = true,
                       renamecols = "_n" => "_np1")

    nrow(joined) == 0 && return nothing

    deltas  = joined.fitness_np1 .- joined.fitness_n
    Δ_abs   = mean(deltas)
    Δ_pct   = 100 * Δ_abs / mean(joined.fitness_n)
    return (Δ_abs, Δ_pct, nrow(joined))
end

# ── Figure ────────────────────────────────────────────────────────────────────

pl = plot(
    layout         = (1, N_CYCLES),
    size           = (340 * N_CYCLES, 500),
    left_margin    = 14mm,
    right_margin   = 6mm,
    top_margin     = 14mm,
    bottom_margin  = 12mm,
    guidefontsize  = 11,
    tickfontsize   = 9,
    titlefontsize  = 12,
    legendfontsize = 8,
)

for cyc in 1:N_CYCLES
    col  = Symbol("n_int_c$(cyc)")
    vals = sort(unique(df[!, col]))

    # Axe Y : clippé au 75e percentile × 2 pour voir la zone utile
    y_max = quantile(df.fitness, 0.75) * 2

    # ── Paired lines : on relie les combinaisons identiques sur les autres cycles
    # Chaque "groupe" = une valeur unique de (n_int_c1, ...) hors cycle courant
    ocols    = other_cols(cyc)
    subdf    = df[!, vcat([col], ocols, [:fitness])]
    grouped  = groupby(subdf, ocols)

    # Palette : une couleur par groupe, triée par fitness à n_interior minimum
    n_groups  = length(grouped)
    palette   = distinguishable_colors(n_groups, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)

    first_group = true
    for (i, g) in enumerate(grouped)
        sg = sort(DataFrame(g), col)
        nrow(sg) < 2 && continue   # pas de ligne si un seul palier couvert

        # Couleur selon la fitness médiane du groupe (pour distinguer bons/mauvais)
        c = palette[i]

        plot!(pl[cyc],
            sg[!, col],
            sg.fitness;
            label             = first_group ? "combinaisons" : "",
            color             = c,
            alpha             = 0.5,
            lw                = 1.5,
            markershape       = :circle,
            markersize        = 5,
            markercolor       = c,
            markerstrokewidth = 0,
        )
        first_group = false
    end

    # ── Courbe du minimum global par palier (en surimpression) ────────────────
    gdf = groupby(df, col)
    agg = combine(gdf, :fitness => minimum => :min_fitness)
    sort!(agg, col)

    plot!(pl[cyc], agg[!, col], agg.min_fitness;
        label             = "min",
        color             = :black,
        lw                = 2.5,
        markershape       = :circle,
        markersize        = 7,
        markercolor       = :black,
        markerstrokewidth = 0,
    )

    # ── Paramètres d'axes ─────────────────────────────────────────────────────
    plot!(pl[cyc];
        ylims  = (0, y_max),
        xlabel = "n_interior  (cycle $cyc)",
        ylabel = cyc == 1 ? "fitness (erreur)" : "",
        title  = "Cycle $cyc",
        xticks = vals,
    )

    # ── Note de lecture ───────────────────────────────────────────────────────
    annotate!(pl[cyc], vals[end], y_max * 0.97,
        text("outliers > $(round(Int, y_max)) masqués", 7, :gray, :right)
    )
end

# ── Sauvegarde ────────────────────────────────────────────────────────────────

mkpath(dirname(OUT_PATH))
savefig(pl, OUT_PATH)
println("Saved → $OUT_PATH")