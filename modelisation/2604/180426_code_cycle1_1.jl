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
using LaTeXStrings


## ===== INPUT DATA =====
df_H = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle1.csv"), DataFrame)
df_V = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle1.csv"), DataFrame)
df_H2 = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle2.csv"), DataFrame)
df_V2 = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle2.csv"), DataFrame)
df_H3 = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle3.csv"), DataFrame)
df_V3 = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle3.csv"), DataFrame)
df_H4 = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle4.csv"), DataFrame)
df_V4 = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle4.csv"), DataFrame)
df_H5 = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle5.csv"), DataFrame)
df_V5 = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle5.csv"), DataFrame)

t_H = df_H[1:6,1]./24; H = df_H[1:6,2]
t_V = df_V[1:6,1]./24; V = df_V[1:6,2]

t_H1 = df_H[1:end,1]./24; H1 = df_H[1:end,2]
t_V1 = df_V[1:end,1]./24; V1 = df_V[1:end,2]

t_H2 = df_H2[1:end,1]./24; H2 = df_H2[1:end,2]
ratio = (t_H2[1]-t_H1[end])
t_H2= t_H2 .- ratio
t_V2 = df_V2[1:end,1]./24; V2 = df_V2[1:end,2]
t_V2= t_V2 .- ratio

t_H3 = df_H3[1:end,1]./24; H3 = df_H3[1:end,2]
ratio = (t_H3[1]-t_H2[end])
t_H3= t_H3 .- ratio
t_V3 = df_V3[1:end,1]./24; V3 = df_V3[1:end,2]
t_V3= t_V3 .- ratio

t_H4 = df_H4[1:end,1]./24; H4 = df_H4[1:end,2]
ratio = (t_H4[1]-t_H3[end])
t_H4= t_H4 .- ratio
t_V4 = df_V4[1:end,1]./24; V4 = df_V4[1:end,2]
t_V4= t_V4 .- ratio

t_H5 = df_H5[1:end,1]./24; H5 = df_H5[1:end,2]
ratio = (t_H5[1]-t_H4[end])
t_H5= t_H5 .- ratio
t_V5 = df_V5[1:end,1]./24; V5 = df_V5[1:end,2]
t_V5= t_V5 .- ratio

t_H_all = vcat(t_H1, t_H2, t_H3, t_H4, t_H5)
H_all   = vcat(H1, H2, H3, H4, H5)

t_V_all = vcat(t_V1, t_V2, t_V3, t_V4, t_V5)
V_all   = vcat(V1, V2, V3, V4, V5)



# CHECK THE DATA ON A PLOTS FIRST
#pl_data = plot(layout=(1,2), size=(700,250), margins=5mm)
#plot!(pl_data[1],t_H,H,label="host data",xlabel="time (days)", ylabel="abundances (cell/ml)", yscale=:log10)
#plot!(pl_data[2],t_V,V,label="virus data",xlabel="time (days)", ylabel="abundances (virus/ml)", yscale=:log10, legend=:bottomright)
# Data looks good we can move on


## ===== CONSTANTS =====
r = 0.5746194091297323
K = 6.675446257207877e7
m = 1.4266424138490254e-11


## ===== MODEL =====
include("model_SV_no_delta.jl")

ϕ = 1.4914407318477944e-8
β = 38.09257390233523

p = [ϕ, β]
u0 = [H[1],V[1]]
tspan = (t_H[1],t_H_all[end])
prob = ODEProblem(model_SV_no_delta, u0, tspan, p)

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
ytick_vals   = [10.0^i for i in 0:2:9]
ytick_labels = [L"10^{%$i}" for i in 0:2:9]

pl_sim = plot(
    layout=(2,1),
    size=(1500,1000),
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
)

# ===== STYLE DATA EXP =====
data_color = RGB(31/255, 119/255, 180/255)
model_color = RGB(255/255, 127/255, 14/255)

# ===== HOST =====
plot!(pl_sim[1], t_H_all, H_all,
    label=" Host data",
    lw=4,
    markershape=:circle,
    markersize=8,
    markerstrokewidth=0,
    color=data_color,
    linealpha=0.7,
    yscale=:log10,
    ylabel="Host abundance\n(cell/mL)"
)

# modèle (ligne continue épaisse)
plot!(pl_sim[1], sol.t, sol[1,:],
    label=" Host model",
    lw=6,
    color=model_color,
    legend=:bottomright
)

# ===== VIRUS =====
plot!(pl_sim[2], t_V_all, V_all,
    label=" Virus data",
    lw=4,
    markershape=:square,
    markersize=8,
    markerstrokewidth=0,
    color=data_color,
    linealpha=0.7,
    yscale=:log10,
    ylabel="Virus abundance\n(particles/mL)"
)

plot!(pl_sim[2], sol.t, sol[2,:],
    label=" Virus model",
    lw=6,
    color=model_color,
    legend=:bottomright
)

# ===== TIMES OF CYCLE CHANGES =====
cycle_changes = [
    t_H2[1],
    t_H3[1],
    t_H4[1],
    t_H5[1]
]

# ===== ADD VERTICAL DOTTED LINES =====
for t_change in cycle_changes
    vline!(pl_sim[1], [t_change],
        color=:grey,
        linestyle=:dot,
        lw=2,
        label=""
    )
    vline!(pl_sim[2], [t_change],
        color=:grey,
        linestyle=:dot,
        lw=2,
        label=""
    )
end

display(pl_sim)