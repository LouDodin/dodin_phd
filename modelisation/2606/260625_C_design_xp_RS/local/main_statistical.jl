using GLMakie
using Printf

# ─────────────────────────────────────────────────────────────────────────────
# PARAMÈTRES
# ─────────────────────────────────────────────────────────────────────────────
const H0      = 1e5
const prop_R0 = 1.0 - 1e-15
const N_R     = prop_R0 * H0   # ≈ 1e5 cellules R par puits


Δt_grid = 10 .^ range(log10(1/24), log10(60.0), length=300)   # 1h → 60j
N_wells_grid = round.(Int, 10 .^ range(log10(1.0), log10(1000.0), length=300))

γ_values = 10 .^ range(log10(1e-5), log10(1e-1), length=300)
γ_ref_idx = argmin(abs.(γ_values .- 9.47815e-3))

# ─────────────────────────────────────────────────────────────────────────────
# FONCTION : calcule P puis la transforme en -log10(1 - P)
# => étire les hautes probabilités (0.9 → 1, 0.99 → 2, 0.999 → 3 …)
# ─────────────────────────────────────────────────────────────────────────────
function compute_logP_matrix(γ)
    M = zeros(length(N_wells_grid), length(Δt_grid))
    for (j, Δt) in enumerate(Δt_grid)
        λ_per_well = γ * N_R * Δt
        for (i, N) in enumerate(N_wells_grid)
            p = 1 - exp(-N * λ_per_well)
            M[i, j] = -log10(1 - p)   # 0 → 0, 0.9 → 1, 0.99 → 2 …
        end
    end
    return M
end

# Ticks de la colorbar : valeurs de P qu'on veut afficher
p_ticks    = [0.0, 0.5, 0.9, 0.95, 0.99, 0.999, 0.9999]
tick_vals  = [-log10(1 - p) for p in p_ticks]
tick_labels = ["0", "50%", "90%", "95%", "99%", "99.9%", "99.99%"]

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE
# ─────────────────────────────────────────────────────────────────────────────
fig = Figure(size = (950, 650))

ax = Axis(fig[1, 1],
    xlabel      = "Δt — durée d'observation (jours)",
    ylabel      = "Nombre de puits",
    title       = "P(observer ≥1 bascule R→S)",
    xscale      = log10,
    yscale      = log10,
    xtickformat = x -> [@sprintf("%.3g", v) for v in x],
    ytickformat = y -> [@sprintf("%d", round(Int, v)) for v in y],
)

# Observable réactif sur γ
γ_obs  = Observable(γ_values[γ_ref_idx])
LP_obs = @lift compute_logP_matrix($γ_obs)

hm = heatmap!(ax,
    Δt_grid, Float64.(N_wells_grid), LP_obs,
    colormap   = :RdYlGn,
    colorrange = (0.0, 4.0),   # 0 → P=0, 4 → P=0.9999
)

Colorbar(fig[1, 2], hm,
    label       = "P(X ≥ 1)",
    ticks       = (tick_vals, tick_labels),
    ticklabelsize = 12,
)

# Contours de référence aux niveaux de P naturels
for (p_level, col, lab) in [
        (0.50, :gray,   "50%"),
        (0.90, :royalblue, "90%"),
        (0.95, :orange, "95%"),
        (0.99, :red,    "99%"),
    ]
    lv = -log10(1 - p_level)
    contour!(ax, Δt_grid, Float64.(N_wells_grid), LP_obs,
        levels    = [lv],
        color     = col,
        linewidth = 2.0,
        linestyle = :dash,
        label     = lab,
    )
end

axislegend(ax, position = :rb, labelsize = 12)

# ── Slider γ ──────────────────────────────────────────────────────────────────
sl_label = Label(fig[2, 1],
    @lift(@sprintf("γ = %.2e j⁻¹", $γ_obs)),
    tellwidth = false, fontsize = 14
)

sl = Slider(fig[3, 1],
    range      = 1:length(γ_values),
    startvalue = γ_ref_idx,
    tellwidth  = true,
)

on(sl.value) do idx
    γ_obs[] = γ_values[idx]
end

rowsize!(fig.layout, 3, Fixed(30))

display(fig)