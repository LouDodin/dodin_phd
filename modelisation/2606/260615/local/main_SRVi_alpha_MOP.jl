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


## ===== Choices =====
const MODEL_FILE = joinpath(@__DIR__, "models/model_SRVi_alpha_MOP.jl")

println("Loading model from: $MODEL_FILE")
include(MODEL_FILE)
using .ModelDef

const prop_S_0 = 1
const V0 = 1e6
const H0 = 5e6
u0 = [prop_S_0 * H0, (1 - prop_S_0) * H0, V0]

model_color = RGB(255/255, 127/255, 14/255)
t = 0:0.1:30
const isoutofdomain = (u, p, t) -> any(x -> x < 0, u)

# ===== Simulation =====

prob = ODEProblem(ModelDef.ODE_MODEL!, u0, (0, 30), nothing)
sol  = solve(prob, Rodas5(),
            reltol=1e-6, abstol=1e-6,
            saveat=t,
            isoutofdomain=isoutofdomain)

# ===== Auto ylims & yticks =====

xtick_vals  = collect(0:3:30)
xtick_labels = string.(collect(0:3:30))

# Subplot 1 : H (= S + R)
H_vals = sol[1, :] .+ sol[2, :]
S_vals = sol[1, :]
R_vals = sol[2, :]
all_vals1 = filter(x -> x > 0, vcat(H_vals, S_vals, R_vals))
exp_min1  = -1
exp_max1  = ceil(Int,  log10(maximum(all_vals1)))
ylims1       = (10.0^exp_min1, 10.0^exp_max1)
ytick_vals1  = [10.0^i for i in exp_min1:exp_max1]
ytick_labels1 = [L"10^{%$i}" for i in exp_min1:exp_max1]

# Subplot 2 : Vi
Vi_vals = sol[3, :]
all_vals2 = filter(x -> x > 0, Vi_vals)
exp_min2  = floor(Int, log10(minimum(all_vals2)))
exp_max2  = ceil(Int,  log10(maximum(all_vals2)))
ylims2       = (10.0^exp_min2, 10.0^exp_max2)
ytick_vals2  = [10.0^i for i in exp_min2:exp_max2]
ytick_labels2 = [L"10^{%$i}" for i in exp_min2:exp_max2]

# ===== Plot =====

pl_fit = plot(
    layout = (2, 1),
    size = (1800, 1000),
    left_margin = 15mm,
    right_margin = 10mm,
    top_margin = 5mm,
    bottom_margin = 10mm,
    grid = true,
    xlims = (0, 30),
    xticks=(xtick_vals, xtick_labels),
    ytickfontsize = 22,
    legendfontsize = 15,
    guidefontsize = 20,
    xtickfontsize = 20,
    titlefontsize = 20,
    xlabel = "Time [days]",
    legend = :bottomright,
    plot_title = "Alpha MOP - H0=$(H0) - V0=$(V0) - prop_S_0=$(prop_S_0)",
    plot_titlefontsize = 25
)

plot!(pl_fit[1], sol.t, H_vals, label="H", color=model_color, lw=4,
      ylabel="Host abundance\n[cells/mL]",
      yscale=:log10, ylims=ylims1,
      yticks=(ytick_vals1, ytick_labels1))
plot!(pl_fit[1], sol.t, S_vals, label="S", color=:red,   lw=2, ls=:dash)
plot!(pl_fit[1], sol.t, R_vals, label="R", color=:green, lw=2, ls=:dash)

plot!(pl_fit[2], sol.t, Vi_vals, label="Vi", color=model_color, lw=4,
      ylabel="Virus abundance\n[parts/mL]",
      yscale=:log10, ylims=ylims2,
      yticks=(ytick_vals2, ytick_labels2))

fmt_sci(x) = begin
    e = floor(Int, log10(x))
    m = x / 10.0^e
    m ≈ 1.0 ? "1e$(e)" : "$(Int(m))e$(e)"
end

display(pl_fit)

mkpath(joinpath(@__DIR__, "output/alpha_MOP"))
fig_path = joinpath(@__DIR__, "output/alpha_MOP/plot_$(fmt_sci(V0))_$(fmt_sci(H0))_$(fmt_sci(prop_S_0)).png")
savefig(pl_fit, fig_path)
println("\nFigure saved to $fig_path")