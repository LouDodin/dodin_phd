## Import packages 
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
using Measures
using LaTeXStrings


## ===== INPUT DATA =====
df_H = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_hostCondition_Temperature20_ReplicateA_Cycle1.csv"), DataFrame)
t_H = df_H[:,1]./24; H = df_H[:,2]

# CHECK THE DATA ON A PLOTS FIRST
pl_data = plot(size=(350,250), margins=5mm)
scatter!(pl_data,t_H,H,label="host data",xlabel="time (days)", ylabel="abundances (cell/ml)", yscale=:log10)
# Data looks good we can move on


## ===== MODEL =====

include("model_S.jl")

r = 0.5746194091297323
K = 6.675446257207877e7
m = 1.4266424138490254e-11

p = [r, K, m]
u0 = [H[1]]
tspan = (t_H[1], t_H[end])
prob = ODEProblem(model_S, u0, tspan, p)

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)

sol = solve(
            prob,
            Tsit5(),
            u0 = u0,
            p = p,
            reltol=1e-6,
            abstol=1e-6,
            isoutofdomain=isoutofdomain
        )

#sol[2,1] = 1E-6

# CHECK THE SIMULATIONS
ytick_vals   = [10.0^i for i in 0:2:8]
ytick_labels = [L"10^{%$i}" for i in 0:2:8]

data_color = RGB(31/255, 119/255, 180/255)
model_color = RGB(255/255, 127/255, 14/255)

pl_sim = plot(
    layout=(2,1),
    size=(1000,1500),
    left_margin=15mm,
    right_margin=10mm,
    top_margin=10mm,
    bottom_margin=10mm,

    grid=true,

    yscale=:log10,
    ylims=(1,2e9),
    yticks=(ytick_vals, ytick_labels),
    ytickfontsize = 28,
    
    legendfontsize=25,
    guidefontsize=25,
    xtickfontsize=25,
    titlefontsize=25,
    
    xlabel="Time (days)",
    ylabel="Host abundance\n(cell/mL)"
)

plot!(pl_sim[1], t_H, H, label=" Host data", markershape=:circle, lw=4, markersize=12, markerstrokewidth=0, linealpha=0.7, color=data_color)
plot!(pl_sim[1], sol.t, sol[1,:], label=" Host model", lw=6, legend=:bottomright, color=model_color)