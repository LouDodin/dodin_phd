## ===== Packages =====
using DifferentialEquations
using OrdinaryDiffEqRosenbrock
using CSV
using DataFrames
using Statistics
using BlackBoxOptim
using Random
using SciMLBase


## ===== Model =====
const MODEL_FILE = joinpath(@__DIR__, "models/SR_RS.jl")
include(MODEL_FILE)
using .ModelDef


## ===== Choices =====
const replicates = ["A", "B", "C"]
const cycles_fit = 5
const cycles_sim = 5
const n_runs     = 1

const lower_bounds = [log(fp.lower) for fp in ModelDef.FITTED_PARAMS]
const upper_bounds = [log(fp.upper) for fp in ModelDef.FITTED_PARAMS]
const isoutofdomain = (u, p, t) -> any(x -> x < 0, u)


## ===== Data =====

raw_data = Dict{String, Vector{NamedTuple}}()

for rep in replicates
    entries    = NamedTuple[]
    t_offset   = nothing
    t_end_prev = nothing

    for cyc_idx in 1:cycles_sim
        path_H = "modelisation/input/xp_input_20/hostData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cyc_idx).csv"
        path_V = "modelisation/input/xp_input_20/virusData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cyc_idx).csv"
        (!isfile(path_H) || !isfile(path_V)) && continue

        df_H = CSV.read(path_H, DataFrame)
        df_V = CSV.read(path_V, DataFrame)

        tH = df_H[:, 1] ./ 24.0;  H = Vector{Float64}(df_H[:, 2])
        tV = df_V[:, 1] ./ 24.0;  V = Vector{Float64}(df_V[:, 2])

        if t_offset === nothing; t_offset = tH[1]; end
        tH .-= t_offset;  tV .-= t_offset

        if t_end_prev !== nothing
            gap = tH[1] - t_end_prev
            tH .-= gap;  tV .-= gap
        end
        t_end_prev = tH[end]

        push!(entries, (index=cyc_idx, tH=tH, H=H, tV=tV, V=V))
    end
    raw_data[rep] = entries
end


## ===== Objective =====
# Each replicate independently (as in SR_RS original)

function objective(θ)
    p = exp.(θ)
    total_err = 0.0

    for rep in replicates
        prop_S = 1.0

        for cycle in 1:cycles_fit
            data = raw_data[rep][cycle]
            u0   = ModelDef.INITIAL_CONDITION(data.H[1], data.V[1], prop_S)

            t_data = sort(unique(vcat(data.tH, data.tV)))
            prob   = ODEProblem(ModelDef.ODE_MODEL!, u0, (t_data[1], t_data[end]), p)
            sol    = solve(prob, Rodas5(), reltol=1e-6, abstol=1e-6,
                           saveat=t_data, isoutofdomain=isoutofdomain)

            if sol.retcode != SciMLBase.ReturnCode.Success || any(x -> x < 0, reduce(vcat, sol.u))
                return 1e15
            end

            H_pred = [sol(t)[1] + sol(t)[2] for t in data.tH]
            V_pred = [sol(t)[3]              for t in data.tV]

            total_err += sum((log.(H_pred) .- log.(data.H)).^2) / length(data.H)
            total_err += sum((log.(V_pred) .- log.(data.V)).^2) / length(data.V)

            u_end  = sol.u[end]
            H_end  = u_end[1] + u_end[2]
            prop_S = H_end > 0 ? clamp(u_end[1] / H_end, 0.0, 1.0) : 1.0
        end
    end

    return total_err
end


## ===== Optimisation =====

results = let
    out = Vector{NamedTuple}(undef, n_runs)
    Threads.@threads for i in 1:n_runs
        Random.seed!(1000 + i)
        res = bboptimize(
            objective;
            SearchRange          = collect(zip(lower_bounds, upper_bounds)),
            NumDimensions        = length(lower_bounds),
            Method               = :adaptive_de_rand_1_bin_radiuslimited,
            PopulationSize       = 250,
            MaxSteps             = 200_000,
            DifferentialWeight   = 0.9,
            crossoverProbability = 0.9,
            TraceMode            = :compact,
        )
        out[i] = (fitness=best_fitness(res), θ=best_candidate(res))
    end
    out
end

p_best = exp.(results[argmin(r.fitness for r in results)].θ)


## ===== Simulate cycle 1 per replicate & export CSV =====

out_dir = joinpath(@__DIR__, "output/$(ModelDef.MODEL_NAME)/cycle1_replicates")
mkpath(out_dir)

for rep in replicates
    data = raw_data[rep][1]
    t0   = min(data.tH[1], data.tV[1])
    t1   = max(data.tH[end], data.tV[end])
    u0   = ModelDef.INITIAL_CONDITION(data.H[1], data.V[1], 1.0)
    t_exp = sort(unique(vcat(data.tH, data.tV)))

    prob = ODEProblem(ModelDef.ODE_MODEL!, u0, (t0, t1), p_best)
    sol  = solve(prob, Rodas5(), reltol=1e-6, abstol=1e-6,
                 saveat=t_exp,
                 isoutofdomain=isoutofdomain)
    
    sol.retcode != SciMLBase.ReturnCode.Success && continue

    S_vec  = [u[1]        for u in sol.u]
    R_vec  = [u[2]        for u in sol.u]
    V_vec  = [u[3]        for u in sol.u]
    SR_vec = S_vec .+ R_vec

    CSV.write(
        joinpath(out_dir, "cycle1_rep$(rep).csv"),
        DataFrame(
            t         = sol.t,
            S         = S_vec,
            R         = R_vec,
            V         = V_vec,
            R_over_SR = [SR > 0 ? R_vec[i] / SR : NaN for (i, SR) in enumerate(SR_vec)]
        )
    )
end