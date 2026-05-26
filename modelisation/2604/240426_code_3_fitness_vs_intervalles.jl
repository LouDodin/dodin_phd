## ===== Import packages =====
using Dates
using CSV
using DataFrames
using Optim
using DifferentialEquations
using LogExpFunctions
using Statistics
using Plots
using Measures
using BlackBoxOptim
using Base.Threads
using DataInterpolations
using Printf
using LaTeXStrings

n_intervalles = [2, 3, 4]
fitness = [36.1, 11.2, 10.1]

## ===== PLOT GLOBAL : fitness vs n_intervalles =====
println("\n===== Génération du plot global fitness vs n_intervalles =====")

pl_global = plot(
    n_intervalles, fitness,
    seriestype = :scatter,
    markershape = :circle,
    markersize = 12,
    markerstrokewidth = 0,
    color = RGB(31/255, 119/255, 180/255),
    xlabel = "Number of intervals per cycle",
    ylabel = "Fitness",
    legend = false,
    grid   = true,
    size   = (1000, 600),
    left_margin  = 12mm,
    right_margin = 8mm,
    top_margin   = 10mm,
    bottom_margin = 10mm,
    xtickfontsize = 18,
    ytickfontsize = 18,
    guidefontsize = 20,
    xticks = n_intervalles,
)

plot!(pl_global, n_intervalles, fitness,
    lw = 2, color = RGB(31/255, 119/255, 180/255), linestyle = :dash, label = nothing)

# Annotations fitness au-dessus de chaque point
for (xi, yi) in zip(n_intervalles, fitness)
    annotate!(pl_global, xi, yi + 2, text(@sprintf("%.3g", yi), 14, :center))
end

global_plot_path = joinpath(@__DIR__, "240426_output/fitness_vs_nintervalles.png")
println("Saving global plot → $global_plot_path")
savefig(pl_global, global_plot_path)

println("\nDone. Summary:")
for (ni, fit) in zip(n_intervalles, fitness)
    @printf("  n_intervalles = %d  →  best_fitness = %.6e\n", ni, fit)
end