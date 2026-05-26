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
using Base.Threads

println(Threads.nthreads())


## ===== INPUT DATA =====
df_H = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle5.csv"), DataFrame)
df_V = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle5.csv"), DataFrame)
t_H = df_H[:,1]./24; H = df_H[:,2]
t_V = df_V[:,1]./24; V = df_V[:,2]

# CHECK THE DATA ON A PLOT FIRST
pl_data = plot(layout=(1,2), size=(700,250), margins=5mm)
scatter!(pl_data[1],t_H,H,label="host data",xlabel="time (days)", ylabel="abundances (cell/ml)", yscale=:log10)
scatter!(pl_data[2],t_V,V,label="virus data",xlabel="time (days)", ylabel="abundances (virus/ml)", yscale=:log10, legend=:bottomright)
# Data looks good we can move on


## ===== CONSTANTS =====
r = 0.5746194091297323
K = 6.675446257207877e7
m = 1.4266424138490254e-11
β = 38.09257390233523
ϕ = 4.28254103565e-10


## ===== SIMULATE =====
include("model_SV_no_delta_only_phi.jl")

p = [ϕ]
u0 = [H[1],V[1]]
tspan = (t_H[1],t_H[end])
prob = ODEProblem(model_SV_no_delta_only_phi, u0, tspan, p)

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)

sol = solve(
            prob,
            Rodas5(),
            u0 = u0,
            p = p,
            reltol=1e-6,
            abstol=1e-6,
            isoutofdomain=isoutofdomain
        )

# CHECK THE SIMULATIONS
pl_sim = plot(layout=(1,2), size=(700,250), margins=5mm)
scatter!(pl_sim[1],t_H,H,label="host data",xlabel="time (days)", ylabel="abundances (cell/ml)", yscale=:log10)
scatter!(pl_sim[2],t_V,V,label="virus data",xlabel="time (days)", ylabel="abundances (virus/ml)", yscale=:log10,legend=:topleft)
plot!(pl_sim[1],sol.t,sol[1,:],label="H",xlabel="time (days)", ylabel="abundances (cell/ml)", yscale=:log10, lw=3, legend=:bottomright)
plot!(pl_sim[2],sol.t,sol[2,:],label="V",xlabel="time (days)", ylabel="abundances (virus/ml)", yscale=:log10, lw=3)