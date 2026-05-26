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

println(Threads.nthreads())


## ===== INPUT DATA =====
df_H = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle5.csv"), DataFrame)
df_V = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle5.csv"), DataFrame)
t_H = df_H[3:end,1]./24; H = df_H[3:end,2]
t_V = df_V[3:end,1]./24; V = df_V[3:end,2]


## ===== CONSTANTS =====
r = 0.5746194091297323
K = 6.675446257207877e7
m = 1.4266424138490254e-11
β = 38.09257390233523


## ===== MODEL =====
include("model_SV_no_delta_only_phi.jl")

p = [1e-10]
u0 = [1e6,1e7]
tspan = (t_H[1],t_H[end])
prob = ODEProblem(model_SV_no_delta_only_phi, u0, tspan, p)

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)


## ===== FITS =====

Y0 = [H[1],V[1]]
tspan = (t_H[1],t_H[end])

function objective_DE(θ)

    # θ is already in log-space
    p = exp.(θ)

    prob_new = remake(prob, u0=Y0, tspan=tspan, p=p)

    sol = solve(
        prob_new,
        Rodas5(),
        reltol=1e-6,
        abstol=1e-6,
        saveat = sort(unique(vcat(t_H, t_V)))
    )

    # reject bad solves
    if sol.retcode != :Success || any(u -> any(x -> !isfinite(x) || x < 0, u), sol.u)
        return 1e12
    end

    YH = sol.(t_H)
    YV = sol.(t_V)

    S_pred = [max(y[1], 1e-12) for y in YH]
    V_pred = [max(y[2], 1e-12) for y in YV]

    err = sum((log.(S_pred) .- log.(H)).^2) +
          sum((log.(V_pred) .- log.(V)).^2)

    return err
end

lower = log.([1E-15])
upper = log.([1E-6])

search_range = [(lower[i], upper[i]) for i in eachindex(lower)]


## VERSION WITH MULTIPLE Run
function run_DE(seed)
    res = bboptimize(
        objective_DE;
        SearchRange = search_range,
        NumDimensions = length(search_range),

        Method = :xnes,

        PopulationSize = 1000,
        MaxSteps = 10000,

        DifferentialWeight = 0.5,
        CrossoverProbability = 0.9,
        
        TraceMode = :silent,   # silent obligatoire en parallèle
        RandomSeed = seed
    )

    return (
        fitness = best_fitness(res),
        θ = best_candidate(res)
    )
end

n_runs = 100

results = Vector{NamedTuple{(:fitness, :θ), Tuple{Float64, Vector{Float64}}}}(undef, n_runs)

Threads.@threads for i in 1:n_runs
    println("Thread ", threadid(), " démarre run ", i)
    results[i] = run_DE(1000 + i)
end

best_idx = argmin(r.fitness for r in results)
best_result = results[best_idx]

θbest = best_result.θ
println("Best error = ", best_result.fitness)
println("Best params = ", exp.(θbest))


#fits
prob_new = remake(prob, u0=Y0, tspan=tspan, p=exp.(θbest))

sol = solve(
        prob_new,
        Rodas5(),
        reltol=1e-6,
        abstol=1e-6,
        isoutofdomain=isoutofdomain
    )


pl_data = plot(layout=(1,2), size=(700,250), margins=5mm)
scatter!(pl_data[1],t_H,H,label="host data",xlabel="time (days)", ylabel="abundances (cell/ml)", yscale=:log10)
scatter!(pl_data[2],t_V,V,label="virus data",xlabel="time (days)", ylabel="abundances (virus/ml)", yscale=:log10,legend=:bottomright)
plot!(pl_data[1],sol.t,sol[1,:],label="host model",lw=3)
plot!(pl_data[2],sol.t,sol[end,:],label="virus model",xlabel="time (days)", ylabel="abundances (virus/ml)",legend=:bottomright,yscale=:log10,lw=3)

println("Saving plot...")
savefig(pl_data, joinpath(@__DIR__, "170426_output/cycle5_2_model_SV_only_phi_no_delta.png"))