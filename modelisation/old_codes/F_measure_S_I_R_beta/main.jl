using DifferentialEquations
using Plots
using Measures
using LaTeXStrings

include(joinpath(@__DIR__, "../models_list.jl"))

host_palette = cgrad([:deepskyblue4, :cyan])
col_S = host_palette[0.2]
col_I = host_palette[0.5]
col_R = host_palette[0.8]
col_H = :deepskyblue4
col_V = :darkred

ti = 7.0
tispan = [0.0, ti]
tf = 40.0
tfspan = [ti, tf]
t = range(0.0, tf, length=Int(tf*100))

modelS = MODELS["S"]
modelI = MODELS["I"]
modelR = MODELS["R"]
model1 = MODELS["SIVi"]
model2 = MODELS["SIRVi_SR"]

const μ   = 0.538001865245161
const k   = 9.784708604680645e7
const φi  = 1.4247469172691298e-8
const β   = 69.99999999999996
const δ   = 0.00010000000056657746
const η   = 5.606949463237242
const εdp = 0.66
const σdp = 0.26
const μ_r = 0.5503400907580093
const k_r = 6.459945866887754e7
const ν   = 1
const α   = 2.4990865344929837e-5

θ_full = (μ, k, φi, β, δ, η, εdp, σdp, μ_r, k_r, ν, α)

θ_S      = (μ, k)
θ_I      = (μ, k, η)
θ_R      = (μ, k, μ_r, k_r)
θ_SIVi   = (μ, k, φi, β, δ, η)
θ_SIRVi  = (μ, k, φi, β, δ, η, μ_r, k_r, α)

floor_plot = 1e-12
floor_plot2 = 1e-5

yticks_vals   = [10.0^i for i in -12:1:9]
yticks_labels = [L"10^{%$i}" for i in -12:1:9]


# =====================
# Well 1: S=1, I=0, R=0
# =====================

X0 = log.([1, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12])

prob1 = ODEProblem(modelS.dynamics!, X0, tispan, θ_S)
sol1 = solve(prob1, Tsit5(), saveat=t[t .<= ti])
X1 = exp.(Array(sol1)')

Xi = copy(X1[end, :])
Xi[4] = (Xi[1] + Xi[2] + Xi[3])*10
Xi = log.(Xi)

prob2 = ODEProblem(model2.dynamics!, Xi, tfspan, θ_SIRVi)
sol2 = solve(prob2, Tsit5(), saveat=t[t .> ti])
X2 = exp.(Array(sol2)')

t_full_1 = vcat(sol1.t, sol2.t)
X_1 = vcat(X1, X2)


# =====================
# Well 2: S=0, I=1, R=0
# =====================

X0 = log.([1e-12, 1, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12])

prob1 = ODEProblem(modelI.dynamics!, X0, tispan, θ_I)
sol1 = solve(prob1, Tsit5(), saveat=t)
X1 = exp.(Array(sol1)')

t_full_2 = sol1.t
X_2 = X1


# =====================
# Well 3: S=0, I=0, R=1
# =====================

X0 = log.([1e-30, 1e-12, 1, 1e-12, 1e-12, 1e-12, 1e-12])

prob1 = ODEProblem(modelR.dynamics!, X0, tispan, θ_R)
sol1 = solve(prob1, Tsit5(), saveat=t[t .<= ti])
X1 = exp.(Array(sol1)')

Xi = copy(X1[end, :])
Xi[4] = (Xi[1] + Xi[2] + Xi[3])*10
Xi = log.(Xi)

prob2 = ODEProblem(model2.dynamics!, Xi, tfspan, θ_SIRVi)
sol2 = solve(prob2, Tsit5(), saveat=t[t .> ti])
X2 = exp.(Array(sol2)')

t_full_3 = vcat(sol1.t, sol2.t)
X_3 = vcat(X1, X2)


# =====================
# Plot
# =====================

plt = plot(layout=(2, 2),
           grid=true,
           yscale=:log10,
           size=(1800, 1200),
           margins=15mm,
           legendfontsize=14,
           guidefontsize=14,
           tickfontsize=14,
           titlefontsize=18)

# --- Panel 1 ---
#plot!(plt[1], t_full_1, max.(X_1[:,1]+X_1[:,2]+X_1[:,3], floor_plot), lw=4, color=col_S, label="H", yticks=(yticks_vals, yticks_labels))
plot!(plt[1], t_full_1, max.(X_1[:,1], floor_plot), lw=4, color=col_S, label="S", yticks=(yticks_vals, yticks_labels))
plot!(plt[1], t_full_1, max.(X_1[:,2], floor_plot), lw=4, color=col_I, label="I")
plot!(plt[1], t_full_1, max.(X_1[:,3], floor_plot), lw=4, color=col_R, label="R")
plot!(plt[1], t_full_1, max.(X_1[:,4] .+ X_1[:,5] .+ X_1[:,6], floor_plot), lw=4, color=col_V, label="Vi")
vline!(plt[1], [ti], lw=2, color=:black, linestyle=:dash, label="ti")
xlabel!(plt[1], "Time (day)")
ylabel!(plt[1], "Concentration (parts/mL)")
title!(plt[1], "Well S: S=1, I=0, R=0")
plot!(plt[1], legend=:bottomright)

# --- Panel 2 ---
plot!(plt[2], t_full_2, max.(X_2[:,1], floor_plot2), lw=4, color=col_S, label="S", yticks=([10.0^i for i in -5:1:2], [L"10^{%$i}" for i in -5:1:2]))
plot!(plt[2], t_full_2, max.(X_2[:,2], floor_plot2), lw=4, color=col_I, label="I")
plot!(plt[2], t_full_2, max.(X_2[:,3], floor_plot2), lw=4, color=col_R, label="R")
plot!(plt[2], t_full_2, max.(X_2[:,4] .+ X_2[:,5] .+ X_2[:,6], floor_plot2), lw=4, color=col_V, label="Vi")
xlabel!(plt[2], "Time (day)")
ylabel!(plt[2], "Concentration (parts/mL)")
title!(plt[2], "Well I: S=0, I=1, R=0")
plot!(plt[2], legend=:bottomright)

# --- Panel 3 ---
plot!(plt[3], t_full_3, max.(X_3[:,1], floor_plot), lw=4, color=col_S, label="S", yticks=(yticks_vals, yticks_labels))
plot!(plt[3], t_full_3, max.(X_3[:,2], floor_plot), lw=4, color=col_I, label="I")
plot!(plt[3], t_full_3, max.(X_3[:,3], floor_plot), lw=4, color=col_R, label="R")
plot!(plt[3], t_full_3, max.(X_3[:,4] .+ X_3[:,5] .+ X_3[:,6], floor_plot), lw=4, color=col_V, label="Vi")
vline!(plt[3], [ti], lw=2, color=:black, linestyle=:dash, label="ti")
xlabel!(plt[3], "Time (day)")
ylabel!(plt[3], "Concentration (parts/mL)")
title!(plt[3], "Well R: S=0, I=0, R=1")
plot!(plt[3], legend=:bottomright)

display(plt)