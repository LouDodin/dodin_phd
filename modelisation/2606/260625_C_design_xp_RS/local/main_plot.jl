using CSV
using DataFrames
using Plots
using LaTeXStrings
using Measures

# lecture du fichier
df = CSV.read(joinpath(@__DIR__, "median_error_vs_gamma.csv"), DataFrame)

γ_theory = df.gamma_true[1]

xtick_vals   = [10.0^i for i in -5:-1]
xtick_labels = [L"10^{%$i}" for i in -5:-1]

# ─────────────────────────────────────────────────────────────────────────────
# SUBPLOT 1 : gamma_median vs gamma_true
# ─────────────────────────────────────────────────────────────────────────────
p1 = scatter(df.gamma_true,
             df.gamma_median_predicted,
             xscale = :log10,
             yscale = :log10,
             marker = :circle,
             xlabel = "γ vrai",
             ylabel = "γ prédit médian",
             title = "γ médian vs γ vrai", label="", legend=:topleft,
             margin=10mm, xticks = (xtick_vals, xtick_labels), yticks = (xtick_vals, xtick_labels))

# ligne y = x
minv = minimum(df.gamma_true)
maxv = maximum(df.gamma_true)

plot!(p1,
      [minv, maxv],
      [minv, maxv],
      lw = 2,
      linestyle = :dash,
      color = :black, alpha=0.5, label = "y=x")

# ligne verticale à γ théorique
vline!(p1, [γ_theory],
       lw = 3,
       linestyle = :dash,
       color = :red, label="γ théorique")
# ─────────────────────────────────────────────────────────────────────────────
# SUBPLOT 2 : erreur médiane
# ─────────────────────────────────────────────────────────────────────────────
p2 = scatter(df.gamma_true,
             df.median_error_percent,
             xscale = :log10,
             marker = :circle,
             xlabel = "γ vrai (log scale)",
             ylabel = "Erreur médiane (%)",
             title = "Erreur médiane vs γ vrai", label="",
             xticks = (xtick_vals, xtick_labels), margin=10mm)

# ligne verticale à γ théorique
vline!(p2, [γ_theory],
       lw = 3,
       linestyle = :dash,
       color = :red, label="γ théorique")

# ─────────────────────────────────────────────────────────────────────────────
# COMBINAISON DES SUBPLOTS
# ─────────────────────────────────────────────────────────────────────────────
p = plot(p1, p2, layout = (1, 2), size = (1000, 400))

savefig(p, joinpath(@__DIR__, "gamma_analysis_subplots.png"))
display(p)