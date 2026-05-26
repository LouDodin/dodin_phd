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


## ===== Infos =====
color_A = RGB(0.6, 0.8, 1.0)
color_B = RGB(31/255, 119/255, 180/255)
color_C = RGB(0.0, 0.3, 0.7) 
model_color = RGB(255/255, 127/255,  14/255)


## ===== Chargement des données =====
function load_replicate(rep::String)
    dfs = []
    for cyc in 1:5
        dfH = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cyc).csv"), DataFrame)
        dfV = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cyc).csv"), DataFrame)
        push!(dfs, (dfH=dfH, dfV=dfV))
    end
    return dfs
end

function build_cycles(dfs)
    cycles = []
    t_offset = 0.0

    for (i, df) in enumerate(dfs)
        tH = df.dfH[:,1] ./ 24
        H  = df.dfH[:,2]
        tV = df.dfV[:,1] ./ 24
        V  = df.dfV[:,2]

        if i == 1
            t_offset = tH[1]
            tH = tH .- t_offset
            tV = tV .- t_offset
        else
            prev = cycles[end]
            gap  = tH[1] / 24 - prev.tH[end]   # gap entre cycles
            tH   = tH .- t_offset .- gap
            tV   = tV .- t_offset .- gap
            t_offset += gap
        end

        push!(cycles, (tH=tH, H=H, tV=tV, V=V, u0=[H[1], V[1]]))
    end
    return cycles
end

dfs_A = load_replicate("A")
dfs_B = load_replicate("B")
dfs_C = load_replicate("C")

cycles_A = build_cycles(dfs_A)
cycles_B = build_cycles(dfs_B)
cycles_C = build_cycles(dfs_C)

all_cycles = vcat(cycles_A, cycles_B, cycles_C)  # 15 cycles au total

# Conditions initiales moyennées par cycle (pour un fit "moyen")
general_cycles = []
for i in 1:5
    cycA = cycles_A[i]
    cycB = cycles_B[i]
    cycC = cycles_C[i]
    tH   = cycA.tH
    tV   = cycA.tV
    H0   = mean([cycA.u0[1], cycB.u0[1], cycC.u0[1]])
    V0   = mean([cycA.u0[2], cycB.u0[2], cycC.u0[2]])
    push!(general_cycles, (tH=tH, tV=tV, u0=[H0, V0]))
end

# Vecteurs globaux pour les plots
t_H_A = vcat([c.tH for c in cycles_A]...);  H_A = vcat([c.H  for c in cycles_A]...)
t_V_A = vcat([c.tV for c in cycles_A]...);  V_A = vcat([c.V  for c in cycles_A]...)
t_H_B = vcat([c.tH for c in cycles_B]...);  H_B = vcat([c.H  for c in cycles_B]...)
t_V_B = vcat([c.tV for c in cycles_B]...);  V_B = vcat([c.V  for c in cycles_B]...)
t_H_C = vcat([c.tH for c in cycles_C]...);  H_C = vcat([c.H  for c in cycles_C]...)
t_V_C = vcat([c.tV for c in cycles_C]...);  V_C = vcat([c.V  for c in cycles_C]...)

# CHECK THE DATA
pl_data = plot(layout=(1,2), size=(900,300), margins=5mm, legend=:topright)
scatter!(pl_data[1], t_H_A, H_A, color=color_A, xlabel="time (days)", ylabel="abundances (cell/ml)", yscale=:log10, ylims=(1e2, 1e8))
scatter!(pl_data[2], t_V_A, V_A, color=color_A, xlabel="time (days)", ylabel="abundances (part/ml)", yscale=:log10, ylims=(1e3, 1e10))
display(pl_data)


## ===== Constants =====
r = 0.5592225270686286
K = 7.29695252684594e7
β = 144
δ = 0.02


## ===== MODEL =====
function model(dY, Y, p, t)
    ϕ  = p[1]
    S  = Y[1]
    Vi = Y[2]
    dY[1] = r*S*(1 - S/K) - ϕ*S*Vi
    dY[2] = β*ϕ*S*Vi - δ*Vi
end

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)

# CHECK THE SIMULATIONS
ϕ = 1e-8
p_test  = [ϕ]
u0_test = [1E6, 1E7]
prob_test = ODEProblem(model, u0_test, (0.0, 10.0), p_test)
sol_test  = solve(prob_test, Tsit5(), reltol=1e-6, abstol=1e-6, isoutofdomain=isoutofdomain)

pl_sim = plot(layout=(1,2), size=(700,250), margins=5mm)
plot!(pl_sim[1], sol_test.t, sol_test[1,:], label="H", xlabel="time (days)",
      ylabel="abundances (cell/ml)", legend=:bottomright, yscale=:log10, lw=3)
plot!(pl_sim[2], sol_test.t, sol_test[2,:], label="V", xlabel="time (days)",
      ylabel="abundances (virus/ml)", legend=:bottomright, yscale=:log10, lw=3)
display(pl_sim)


## ===== OBJECTIVE (tous les réplicats, tous les cycles) =====
function objective(θ)
    total_err = 0.0
    p_model   = exp.(θ)

    for cyc in all_cycles          # 15 cycles (3 réplicats × 5 cycles)
        t0     = cyc.tH[1]
        t1     = max(cyc.tH[end], cyc.tV[end])
        t_save = sort(unique(vcat(cyc.tH, cyc.tV)))

        prob_cyc = ODEProblem(model, cyc.u0, (t0, t1), p_model)
        sol_cyc  = solve(
            prob_cyc, Rodas5(),
            reltol        = 1e-6,
            abstol        = 1e-6,
            saveat        = t_save,
            isoutofdomain = isoutofdomain,
        )

        if sol_cyc.retcode != ReturnCode.Success ||
           any(u -> any(x -> !isfinite(x) || x < 0, u), sol_cyc.u)
            return 1e12
        end

        S_pred = [max(sol_cyc(t)[1], 1e-12) for t in cyc.tH]
        V_pred = [max(sol_cyc(t)[2], 1e-12) for t in cyc.tV]

        total_err += sum((log.(S_pred) .- log.(cyc.H)).^2) / length(cyc.tH)
        total_err += sum((log.(V_pred) .- log.(cyc.V)).^2) / length(cyc.tV)
    end

    return total_err
end


## ===== MULTI-RUN (parallelised with Threads) =====
lower        = log.([1E-15])
upper        = log.([1E-6])
search_range = [(lower[i], upper[i]) for i in eachindex(lower)]

function run_DE(seed)
    res = Base.invokelatest(
        bboptimize, objective;
        SearchRange          = search_range,
        NumDimensions        = length(search_range),
        Method               = :xnes,
        PopulationSize       = 1000,
        MaxSteps             = 10000,
        DifferentialWeight   = 0.5,
        CrossoverProbability = 0.9,
        TraceMode            = :silent,
        RandomSeed           = seed,
    )
    return (fitness=best_fitness(res), θ=best_candidate(res))
end

n_runs  = 1
results = Vector{NamedTuple}(undef, n_runs)

Threads.@threads for i in 1:n_runs
    println("  Thread ", threadid(), " démarre run ", i, " (seed=$(1000+i))")
    results[i] = Base.invokelatest(run_DE, 1000 + i)
end

best_idx    = argmin(r.fitness for r in results)
best_result = results[best_idx]
θbest       = best_result.θ
println("Best error  = ", best_result.fitness)
println("Best ϕ      = ", exp.(θbest))


## ===== PLOT FINAL =====
# Un solve par cycle avec les conditions initiales de chaque réplicat
prob_ref = ODEProblem(model, all_cycles[1].u0, (0.0, 1.0), exp.(θbest))

pl_fit = plot(layout=(1,2), size=(1200,400), margins=10mm, legend=:topright)

rep_colors = [data_color_A, data_color_B, data_color_C]
rep_labels = ["Rep A", "Rep B", "Rep C"]
rep_cycles = [cycles_A, cycles_B, cycles_C]

for (ri, (cycs, col, lab)) in enumerate(zip(rep_cycles, rep_colors, rep_labels))
    for (ci, cyc) in enumerate(cycs)
        t0   = cyc.tH[1]
        t1   = max(cyc.tH[end], cyc.tV[end])
        sol_c = solve(
            remake(prob_ref, u0=cyc.u0, tspan=(t0, t1), p=exp.(θbest)),
            Rodas5(), reltol=1e-6, abstol=1e-6, isoutofdomain=isoutofdomain
        )
        lbl = ci == 1 ? lab : ""   # légende une seule fois par réplicat
        scatter!(pl_fit[1], cyc.tH, cyc.H,   label=lbl, color=col, alpha=0.6,
                 xlabel="time (days)", ylabel="abundances (cell/ml)",
                 yscale=:log10, ylims=(1e2, 1e8))
        scatter!(pl_fit[2], cyc.tV, cyc.V,   label=lbl, color=col, alpha=0.6,
                 xlabel="time (days)", ylabel="abundances (part/ml)",
                 yscale=:log10, ylims=(1e3, 1e10))
        plot!(pl_fit[1], sol_c.t, sol_c[1,:], label="", color=:black, lw=1.5, alpha=0.5)
        plot!(pl_fit[2], sol_c.t, sol_c[2,:], label="", color=:black, lw=1.5, alpha=0.5)
    end
end

display(pl_fit)
savefig(pl_fit, joinpath(@__DIR__, "290426_output/SV_model.png"))