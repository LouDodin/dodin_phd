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


## ===== INPUT DATA =====
df_A = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_hostCondition_Temperature20_ReplicateA_Cycle1.csv"), DataFrame)
df_B = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_hostCondition_Temperature20_ReplicateB_Cycle1.csv"), DataFrame)
df_C = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_hostCondition_Temperature20_ReplicateC_Cycle1.csv"), DataFrame)

t_A = df_A[:,1]./24; H_A = df_A[:,2]
t_B = df_B[:,1]./24; H_B = df_B[:,2]
t_C = df_C[:,1]./24; H_C = df_C[:,2]

replicates = [
    (t=t_A, H=H_A),
    (t=t_B, H=H_B),
    (t=t_C, H=H_C),
]

t_min   = minimum([rep.t[1]   for rep in replicates])
t_max   = maximum([rep.t[end] for rep in replicates])
H0_mean = mean([rep.H[1]      for rep in replicates])

Y0    = [H0_mean]
tspan = (t_min, t_max)

# CHECK THE DATA
pl_data = plot(size=(350,250), margins=5mm)
scatter!(pl_data, t_A, H_A, label="Rep A", xlabel="time (days)", ylabel="abundances (cell/ml)", yscale=:log10)
scatter!(pl_data, t_B, H_B, label="Rep B")
scatter!(pl_data, t_C, H_C, label="Rep C")
display(pl_data)


## ===== MODEL =====
include("model_S_2.jl")

r = 0.5881765172
K = 6.0e7

p = [r, K]
prob = ODEProblem(model, Y0, tspan, p)

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)


## ===== JOINT OBJECTIVE =====
function objective_joint(θ)
    p = exp.(θ)

    all_t = sort(unique(vcat([rep.t for rep in replicates]...)))

    prob_new = remake(prob, u0=Y0, tspan=tspan, p=p)

    sol = solve(
        prob_new,
        Rodas5(),
        reltol=1e-6,
        abstol=1e-6,
        saveat=all_t
    )

    if sol.retcode != :Success || any(u -> any(x -> !isfinite(x) || x < 0, u), sol.u)
        return 1e12
    end

    total_err = 0.0
    for rep in replicates
        S_pred = [max(sol(t)[1], 1e-12) for t in rep.t]
        total_err += sum((log.(S_pred) .- log.(rep.H)).^2)
    end

    return total_err
end


## ===== MULTI-RUN ON 100 SEEDS (parallelised with Threads) =====
lower = log.([0.1, 1E7])
upper = log.([0.7, 1E8])
search_range = [(lower[i], upper[i]) for i in eachindex(lower)]

function run_DE(seed)
    res = bboptimize(
        objective_joint;
        SearchRange = search_range,
        NumDimensions = length(search_range),
        Method = :xnes,
        PopulationSize = 1000,
        MaxSteps = 10000,
        DifferentialWeight = 0.5,
        CrossoverProbability = 0.9,
        TraceMode = :silent,
        RandomSeed = seed
    )
    return (
        fitness = best_fitness(res),
        θ = best_candidate(res)
    )
end

n_runs = 100
results = Vector{NamedTuple}(undef, n_runs)

Threads.@threads for i in 1:n_runs
    println("  Thread ", threadid(), " démarre run ", i, " (seed=$(1000+i))")
    results[i] = run_DE(1000 + i)
end

best_idx = argmin(r.fitness for r in results)
best_result = results[best_idx]

θbest = best_result.θ
println("Best error = ", best_result.fitness)
println("Best params = ", exp.(θbest))


## ===== PLOT =====
prob_best = remake(prob, u0=Y0, tspan=tspan, p=exp.(θbest))

sol = solve(
    prob_best,
    Rodas5(),
    reltol=1e-6,
    abstol=1e-6,
    isoutofdomain=isoutofdomain
)

labels = ["Rep A", "Rep B", "Rep C"]
colors = [:blue, :red, :green]

pl_fit = plot(size=(400,300), margins=5mm, xlabel="time (days)", ylabel="abundances (cell/ml)", yscale=:log10)

for (i, rep) in enumerate(replicates)
    scatter!(pl_fit, rep.t, rep.H, label=labels[i], color=colors[i])
end

plot!(pl_fit, sol.t, sol[1,:], label="model", color=:black, lw=2)

display(pl_fit)
savefig(pl_fit, joinpath(@__DIR__, "230426_output/fit_r_K_3rep.png"))