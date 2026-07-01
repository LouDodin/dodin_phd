using GLMakie
using Printf

# ─────────────────────────────────────────────────────────────────────────────
# PARAMÈTRES
# ─────────────────────────────────────────────────────────────────────────────
const H0      = 1e5
const prop_R0 = 1.0 - 1e-15
const N_R     = prop_R0 * H0   # ≈ 1e5 cellules R par puits (population INITIALE)

# Dynamique de la population R :  dR/dt = r*R*(1 - R/K) - γ*R
r = 0.574619342477644
K = 6.675449070379925e7
const R0 = N_R        # condition initiale R(0)

Δt_grid = 10 .^ range(log10(1/24), log10(60.0), length=300)   # 1h → 60j
N_wells_grid = round.(Int, 10 .^ range(log10(1.0), log10(1000.0), length=300))

γ_values = 10 .^ range(log10(1e-5), log10(1e-1), length=300)
γ_ref_idx = argmin(abs.(γ_values .- 9.47815e-3))

# ─────────────────────────────────────────────────────────────────────────────
# DYNAMIQUE R(t) : solution fermée de la logistique avec perte -γR
#   dR/dt = rR(1-R/K) - γR  =  ρR(1-R/K_eff),  ρ = r-γ,  K_eff = K*ρ/r
# ─────────────────────────────────────────────────────────────────────────────
function R_of_t(t::AbstractVector, γ)
    ρ = r - γ
    if abs(ρ) < 1e-10
        # cas dégénéré r ≈ γ : dR/dt = -(r/K) R²
        return R0 ./ (1 .+ (r/K) .* R0 .* t)
    else
        Keff = K * ρ / r
        A = (Keff - R0) / R0
        return Keff ./ (1 .+ A .* exp.(-ρ .* t))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# INTÉGRALE CUMULÉE Λ(Δt) = ∫₀^Δt γ R(t) dt   (processus de Poisson non-homogène)
# calculée par trapèzes sur une grille fine, puis interpolée sur Δt_grid
# ─────────────────────────────────────────────────────────────────────────────
function cumtrapz(x, y)
    n = length(x)
    out = zeros(n)
    @inbounds for i in 2:n
        out[i] = out[i-1] + (x[i] - x[i-1]) * (y[i] + y[i-1]) / 2
    end
    return out
end

function interp_linear(xs, ys, x)
    x <= xs[1]   && return ys[1]
    x >= xs[end] && return ys[end]
    i = searchsortedfirst(xs, x)
    x0, x1 = xs[i-1], xs[i]
    y0, y1 = ys[i-1], ys[i]
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
end

function compute_Lambda_grid(γ; n_fine = 4000)
    t_fine = collect(range(0.0, maximum(Δt_grid), length = n_fine))
    Rt = R_of_t(t_fine, γ)
    integrand = γ .* Rt                      # taux instantané de bascule
    Λ_fine = cumtrapz(t_fine, integrand)      # intégrale cumulée
    return [interp_linear(t_fine, Λ_fine, Δt) for Δt in Δt_grid]
end

# ─────────────────────────────────────────────────────────────────────────────
# FONCTION : calcule P puis la transforme en -log10(1 - P)
# ─────────────────────────────────────────────────────────────────────────────
function compute_logP_matrix(γ)
    Λ_grid = compute_Lambda_grid(γ)   # Λ(Δt) pour chaque colonne
    M = zeros(length(N_wells_grid), length(Δt_grid))
    for (j, Λ) in enumerate(Λ_grid)
        for (i, N) in enumerate(N_wells_grid)
            p = 1 - exp(-N * Λ)
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
    title       = "P(observer ≥1 bascule R→S)  —  dR/dt = rR(1-R/K) - γR",
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
    @lift(@sprintf("γ = %.2e j⁻¹   (r = %.2g j⁻¹, K = %.2g)", $γ_obs, r, K)),
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