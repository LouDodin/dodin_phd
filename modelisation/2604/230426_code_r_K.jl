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


## ===== INPUT DATA =====
df_H = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_hostCondition_Temperature20_ReplicateA_Cycle1.csv"), DataFrame)
t_H = df_H[:,1]./24; H = df_H[:,2]

# CHECK THE DATA ON A PLOTS FIRST
pl_data = plot(size=(350,250), margins=5mm)
scatter!(pl_data,t_H,H,label="host data",xlabel="time (days)", ylabel="abundances (cell/ml)", yscale=:log10)
# Data looks good we can move on


## ===== MODEL =====

include("model_S_2.jl")

r = 0.5881765172
K = 6.0e7

p = [r, K]
u0 = [1E6]
tspan = (0, 10)
prob = ODEProblem(model, u0, tspan, p)

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
pl_sim = plot(size=(350,250), margins=5mm)
plot!(pl_sim,sol.t,sol[1,:],label="H",xlabel="time (days)", ylabel="abundances (cell/ml)", yscale=:log10, lw=3, legend=:bottomleft)


## ===== FITS =====

Y0 = [H[1]]
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
        saveat = t_H
    )

    # reject bad solves
    if sol.retcode != :Success || any(u -> any(x -> !isfinite(x) || x < 0, u), sol.u)
        return 1e12
    end

    YH = sol.(t_H)

    S_pred = [max(y[1], 1e-12) for y in YH]

    err = sum((log.(S_pred) .- log.(H)).^2)

    return err
end

lower = log.([0.1, 1E7])
upper = log.([0.7, 1E8])

search_range = [(lower[i], upper[i]) for i in eachindex(lower)]

res = bboptimize(
    objective_DE;
    SearchRange = search_range,
    NumDimensions = length(search_range),

    Method = :xnes,

    PopulationSize = 1000,
    MaxSteps = 10000,

    DifferentialWeight = 0.5,
    CrossoverProbability = 0.9,

    TraceMode = :verbose,
    RandomSeed = 123
)

θbest = best_candidate(res)
println(exp.(best_candidate(res)))

prob_new = remake(prob, u0=Y0, tspan=tspan, p=exp.(θbest))

sol = solve(
        prob_new,
        Rodas5(),
        reltol=1e-6,
        abstol=1e-6,
        isoutofdomain=isoutofdomain
    )

pl_data = plot(size=(350,250), margins=5mm)
scatter!(pl_data, t_H, H, label="host data", xlabel="time (days)", ylabel="abundances (cell/ml)", yscale=:log10)
plot!(pl_data, sol.t, sol[1,:], label="Hmod", lw=3)


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

        TraceMode = :verbose,
        RandomSeed = seed
    )

    return (
        fitness = best_fitness(res),
        θ = best_candidate(res)
    )
end

n_runs = 100

results = [run_DE(1000 + i) for i in 1:n_runs]

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


pl_data = plot(size=(350,250), margins=5mm)
scatter!(pl_data,t_H,H,label="host data",xlabel="time (days)", ylabel="abundances (cell/ml)", yscale=:log10)
plot!(pl_data,sol.t,sol[1,:],label="host model",lw=3)


println("Saving plot...")
savefig(pl_data, joinpath(@__DIR__, "230426_output/fit_r_K.png"))