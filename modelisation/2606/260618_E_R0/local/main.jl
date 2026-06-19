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
using Colors


## ===== Choices =====
const MODEL_FILE = joinpath(@__DIR__, "models/SR_RS.jl")
const replicates  = ["A", "B", "C"]
const cycles_fit  = 5
const cycles_sim  = 5
const n_runs      = 1


## ===== Utils =====
include(MODEL_FILE)
using .ModelDef

lower_bounds = [log(fp.lower) for fp in ModelDef.FITTED_PARAMS]
upper_bounds = [log(fp.upper) for fp in ModelDef.FITTED_PARAMS]

output_dir = joinpath(@__DIR__, "output/$(ModelDef.MODEL_NAME)")
mkpath(output_dir)

const isoutofdomain = (u, p, t) -> any(x -> x < 0, u)

replicate_colors = Dict(
    "A" => RGB(0.6, 0.8, 1.0),
    "B" => RGB(31/255, 119/255, 180/255),
    "C" => RGB(0.0, 0.3, 0.7)
)

model_color = RGB(255/255, 127/255, 14/255)
cycle_color = RGB(31/255, 119/255, 180/255)


## ===== Input: S & V data =====

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

# Check plot
pl_data = plot(layout=(1,2), size=(900,350), margins=5mm, legend=:bottomright)

for rep in replicates
    for cyc in 1:cycles_sim
        data = raw_data[rep][cyc]
        scatter!(pl_data[1], data.tH, data.H;
            color=replicate_colors[rep], label=cyc==1 ? "Rep $rep" : "",
            xlabel="Time (days)", ylabel="Host abundance (cell/mL)",
            yscale=:log10, ylims=(1e2,1e8), title="H")
        scatter!(pl_data[2], data.tV, data.V;
            color=replicate_colors[rep], label=cyc==1 ? "Rep $rep" : "",
            xlabel="Time (days)", ylabel="Virus abundance (part/mL)",
            yscale=:log10, ylims=(1e3,1e10), title="Vi")
    end
end
#display(pl_data)


## ===== Helper functions =====

# --- ODE solver ---
function solve_cycle(p, u0, t0, t1, t_save)
    prob = ODEProblem(ModelDef.ODE_MODEL!, u0, (t0, t1), p)
    sol  = solve(prob,
                Rodas5(),
                reltol=1e-6,
                abstol=1e-6,
                saveat=t_save,
                isoutofdomain=isoutofdomain
        )
    return sol
end

# --- Objective function ---
# Each replicate independently

function objective(θ)
    p = exp.(θ)
    prop_S0 = p[2]

    total_err = 0.0

    for rep in replicates

        prop_S  = prop_S0  # fraction of susceptible hosts at cycle start

        for cycle in 1:cycles_fit
            data = raw_data[rep][cycle]
            H0 = data.H[1]
            V0 = data.V[1]
            u0 = ModelDef.INITIAL_CONDITION(H0, V0, prop_S)

            t_data = sort(unique(vcat(data.tH, data.tV)))
            t0_data, t1_data = t_data[1], t_data[end]

            sol = solve_cycle(p, u0, t0_data, t1_data, t_data)

            if sol.retcode != SciMLBase.ReturnCode.Success || any(x -> x < 0, reduce(vcat, sol.u))
                return 1e15
            end

            H_pred = [sol(t)[1]+sol(t)[2] for t in data.tH]
            V_pred = [sol(t)[3] for t in data.tV]

            total_err += sum((log.(H_pred) .- log.(data.H)).^2) / length(data.H)
            total_err += sum((log.(V_pred) .- log.(data.V)).^2) / length(data.V)
            
            # Proportions for next cycle
            u_end   = sol.u[end]
            H_end   = u_end[1] + u_end[2]
            prop_S  = u_end[1] / H_end
        end
    end

    return total_err
end


## ===== Optimisation =====
function run_DE(seed)
    Random.seed!(seed)
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
    return (fitness=best_fitness(res), θ=best_candidate(res))
end

function main()
    results = Vector{NamedTuple}(undef, n_runs)
    Threads.@threads for i in 1:n_runs
        println("  Thread $(Threads.threadid()) — run $i (seed=$(1000+i))")
        results[i] = run_DE(1000 + i)
    end
    return results
end

results = Base.invokelatest(main)

best_idx    = argmin(r.fitness for r in results)
best_result = results[best_idx]
θbest       = best_result.θ
p_best      = exp.(θbest)

println("\n===== Best result ($(ModelDef.MODEL_NAME)) =====")
println("  Fitness = ", best_result.fitness)
for (i, fp) in enumerate(ModelDef.FITTED_PARAMS)
    println("  $(fp.name) = $(p_best[i])  ($(fp.description))")
end
println("  All run fitnesses: ", [r.fitness for r in results])


## ===== Information Criteria =====

# --- 1. Count total observations ---
n_obs = sum(
    length(raw_data[rep][cyc].H) + length(raw_data[rep][cyc].V)
    for rep in replicates
    for cyc in 1:cycles_fit
)

# --- 2. Recompute raw RSS (sum, not mean) from best parameters ---
function raw_rss(θ)
    p = exp.(θ)
    prop_S0 = p[2] 
    rss = 0.0

    for rep in replicates
        prop_S = prop_S0
        for cycle in 1:cycles_fit
            data = raw_data[rep][cycle]
            H0   = data.H[1];  V0 = data.V[1]
            u0   = ModelDef.INITIAL_CONDITION(H0, V0, prop_S)

            t_data = sort(unique(vcat(data.tH, data.tV)))
            sol    = solve_cycle(p, u0, t_data[1], t_data[end], t_data)

            (sol.retcode != SciMLBase.ReturnCode.Success ||
             any(x -> x < 0, reduce(vcat, sol.u))) && return Inf

            H_pred = [sol(t)[1] + sol(t)[2] for t in data.tH]
            V_pred = [sol(t)[3] for t in data.tV]

            rss += sum((log.(H_pred) .- log.(data.H)).^2)
            rss += sum((log.(V_pred) .- log.(data.V)).^2)

            u_end  = sol.u[end]
            H_end  = u_end[1] + u_end[2]
            prop_S = H_end > 0 ? clamp(u_end[1] / H_end, 0.0, 1.0) : 1.0
        end
    end
    return rss
end

rss_best = raw_rss(θbest)

# --- 3. MLE log-likelihood (Gaussian errors on log-scale, σ² estimated as RSS/n) ---
k    = length(ModelDef.FITTED_PARAMS)
σ²   = rss_best / n_obs
logL = -n_obs/2 * log(2π) - n_obs/2 * log(σ²) - n_obs/2

# --- 4. Criteria ---
AIC  = 2k - 2logL
AICc = AIC + 2k*(k+1) / (n_obs - k - 1)   # correction for small n
BIC  = k*log(n_obs) - 2logL

println("\n===== Information Criteria ($(ModelDef.MODEL_NAME)) =====")
println("  n_obs = $n_obs,  k = $k")
println("  RSS   = $(round(rss_best,  sigdigits=4))")
println("  logL  = $(round(logL,      sigdigits=6))")
println("  AIC   = $(round(AIC,       sigdigits=6))")
println("  AICc  = $(round(AICc,      sigdigits=6))")
println("  BIC   = $(round(BIC,       sigdigits=6))")


## ===== Phi_ref =====
const poly_file = "modelisation/2606/260604/genotoul/output/nint_3_2_2_3_2/polynomial.txt"


lines = readlines(poly_file)

function extract_float(s::AbstractString)::Float64
    m = match(r"([+-]?\s*[0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)", strip(s))
    m === nothing && error("Cannot extract float from: \"$s\"")
    parse(Float64, replace(m.captures[1], " " => ""))
end

ϕ_intervals = Vector{NTuple{6,Float64}}()
global i = 1
while i <= length(lines)
    m_iv = match(r"Interval\s+\[([0-9eE+\-.]+),\s*([0-9eE+\-.]+)\]\s+days:", strip(lines[i]))
    if m_iv !== nothing
        if i + 4 > length(lines)
            global i += 1; continue
        end
        t0 = parse(Float64, m_iv.captures[1])
        t1 = parse(Float64, m_iv.captures[2])
        a  = extract_float(split(lines[i+1], "=")[end])
        b  = extract_float(split(lines[i+2], "*")[1])
        c  = extract_float(split(lines[i+3], "*")[1])
        d  = extract_float(split(lines[i+4], "*")[1])
        push!(ϕ_intervals, (t0, t1, a, b, c, d))
        i += 5; continue
    end
    global i += 1
end

const tmin   = ϕ_intervals[1][1]
const tmax   = ϕ_intervals[end][2]

function ϕ_ref(t)
    tc  = clamp(t, tmin, tmax)
    idx = length(ϕ_intervals)
    for k in eachindex(ϕ_intervals)
        if tc <= ϕ_intervals[k][2]; idx = k; break; end
    end
    t0, _, a, b, c, d = ϕ_intervals[idx]
    dt = tc - t0
    exp(a + b*dt + c*dt^2 + d*dt^3)
end

const N_GRID  = 10_000
const t_grid  = collect(range(tmin, tmax; length=N_GRID))
const ϕ_grid = ϕ_ref.(t_grid)

knots = Float64[]
for interval in ϕ_intervals
    push!(knots, interval[1])
end
push!(knots, ϕ_intervals[end][2])
unique!(sort!(knots))


## ===== Simulation & Plot =====
ytick_vals1   = [10.0^i for i in 2:9]
ytick_labels1 = [L"10^{%$i}" for i in 2:9]
ytick_vals2   = [10.0^i for i in 5:10]
ytick_labels2 = [L"10^{%$i}" for i in 5:10]
ytick_vals3   = [10.0^i for i in -14:-7]
ytick_labels3 = [L"10^{%$i}" for i in -14:-7]


param_str = join(
    ["$(ModelDef.FITTED_PARAMS[i].name)=$(round(p_best[i], sigdigits=2))" 
     for i in eachindex(ModelDef.FITTED_PARAMS)],
    "  "
)

pl_fit = plot(
    layout = (3, 1),
    size = (1800, 2000),
    left_margin = 15mm,
    right_margin = 10mm,
    top_margin = 35mm,
    bottom_margin = 5mm,
    grid = true,
    yscale = :log10,
    xlims = (0, 67),
    ytickfontsize = 22,
    legendfontsize = 15,
    guidefontsize = 20,
    xtickfontsize = 20,
    titlefontsize = 20,
    xlabel = "Time [days]",
    legend = :bottomright,
    plot_title = "$(ModelDef.MODEL_NAME)\n$(param_str)\nFitness = $(round(best_result.fitness, sigdigits=4))\nAIC = $(round(AIC, sigdigits=4)) - AICc = $(round(AICc, sigdigits=4)) - BIC = $(round(BIC, sigdigits=4))",
    plot_titlefontsize = 25
)

# --- Data ---
for rep in replicates
    for cyc_idx in 1:cycles_sim
        data = raw_data[rep][cyc_idx]
        lbl  = cyc_idx == 1 ? "Replicate $rep" : ""
        scatter!(pl_fit[1], data.tH, data.H;
            label=lbl, color=replicate_colors[rep], alpha=0.7,
            ylabel="Host\n[cell/mL]",
            ylims=(1e2, 1e9), yticks=(ytick_vals1, ytick_labels1),
            markershape=:circle, markersize=10,
            legend=(0.15, 1))
        scatter!(pl_fit[2], data.tV, data.V;
            label=lbl, color=replicate_colors[rep], alpha=0.7,
            ylabel="Virus\n[part/mL]",
            ylims=(1e5, 1e10), yticks=(ytick_vals2, ytick_labels2),
            markershape=:circle, markersize=10,
            legend=(0.15, 0.3))
    end
end

plot!(pl_fit[3], t_grid, ϕ_grid;
      label="ϕ_ref", color=replicate_colors["B"], lw=6, alpha=0.7,
      ylabel="ϕ\n[mL/(part.day)]",
      ylims=(1e-14, 1e-7), yticks=(ytick_vals3, ytick_labels3))

phi_knot_values = ϕ_ref.(knots)
scatter!(pl_fit[3], knots, phi_knot_values, color=replicate_colors["C"],
      markersize=10, markershape=:diamond,
      label="Knots", legend=:bottomright)


# --- Model ---
global prop_S_plot = p_best[2]

for cycle in 1:cycles_sim
    H0 = mean(raw_data[rep][cycle].H[1] for rep in replicates)
    V0 = mean(raw_data[rep][cycle].V[1] for rep in replicates)
    u0 = ModelDef.INITIAL_CONDITION(H0, V0, prop_S_plot)

    t0 = minimum(raw_data[rep][cycle].tH[1]   for rep in replicates)
    t1 = maximum(max(raw_data[rep][cycle].tH[end], raw_data[rep][cycle].tV[end]) for rep in replicates)

    sol = solve(ODEProblem(ModelDef.ODE_MODEL!, u0, (t0, t1), p_best),
                Rodas5(), reltol=1e-6, abstol=1e-6,
                isoutofdomain=isoutofdomain)

    # Proportions for next cycle
    u_end = sol.u[end]
    H_end = u_end[1] + u_end[2]
    global prop_S_plot  = H_end > 0 ? clamp(u_end[1]  / H_end, 0.0, 1.0) : 1.0
    
    S = [u[1] for u in sol.u]
    R = [u[2] for u in sol.u]
    V = [u[3] for u in sol.u]
    H = S .+ R
    ϕ_equiv = [ModelDef.PHI_EQUIV(sol(t), p_best) for t in sol.t]

    plot!(pl_fit[1], sol.t, H, label= cycle==1 ? "H model" : "", color=model_color, lw=6, alpha=0.7)
    plot!(pl_fit[1], sol.t, S, label= cycle==1 ? "S model" : "", color=:red, lw=4, ls=:dash, alpha=0.7)
    plot!(pl_fit[1], sol.t, R, label= cycle==1 ? "R model" : "", color=:green, lw=4, ls=:dash, alpha=0.7)
    plot!(pl_fit[2], sol.t, V, label= cycle==1 ? "V model" : "", color=model_color, lw=6, alpha=0.7)
    plot!(pl_fit[3], sol.t, ϕ_equiv, label= cycle==1 ? "ϕ_equiv model" : "", color=model_color, lw=6, alpha=0.7)
end

# --- Dilutions ---
cycle_tbounds = Vector{Tuple{Float64,Float64}}(undef, cycles_sim)
for cycle in 1:cycles_sim
    t0 = minimum(raw_data[rep][cycle].tH[1] for rep in replicates)
    t1 = maximum(max(raw_data[rep][cycle].tH[end],
                     raw_data[rep][cycle].tV[end]) for rep in replicates)
    cycle_tbounds[cycle] = (t0, t1)
end

dilution_times = [cycle_tbounds[c][2] for c in 1:cycles_sim-1]

for panel in 1:3
    for (i, td) in enumerate(dilution_times)
        vline!(pl_fit[panel], [td];
               color=cycle_color, lw=2, ls=:dot,
               label = nothing)
    end
end

display(pl_fit)

fig_path = "$(output_dir)/$(ModelDef.MODEL_NAME)_plot.png"
savefig(pl_fit, fig_path)
println("\nFigure saved to $fig_path")


## ===== Log =====
log_path = "$(output_dir)/$(ModelDef.MODEL_NAME)_log.txt"
open(log_path, "w") do io
    println(io, "===== Run log =====")
    println(io, "Date        : ", Dates.now())
    println(io, "Model       : ", ModelDef.MODEL_NAME)
    println(io, "Model file  : ", MODEL_FILE)
    println(io, "n_runs      : ", n_runs)
    println(io, "cycles_fit  : ", cycles_fit)
    println(io, "cycles_sim  : ", cycles_sim)
    println(io, "replicates  : ", replicates)
    println(io)
    println(io, "===== Fixed parameters =====")
    for (k, v) in pairs(ModelDef.FIXED_PARAMS)
        println(io, "  $k = $v")
    end
    println(io)
    println(io, "===== Search bounds =====")
    for fp in ModelDef.FITTED_PARAMS
        println(io, "  $(fp.name)  ∈ [$(fp.lower), $(fp.upper)]  — $(fp.description)")
    end
    println(io)
    println(io, "===== All runs =====")
    for (i, res) in enumerate(results)
        p_i = exp.(res.θ)
        param_str = join(["$(ModelDef.FITTED_PARAMS[j].name)=$(p_i[j])" for j in eachindex(ModelDef.FITTED_PARAMS)], "  ")
        println(io, "  Run $i (seed=$(1000+i)) : loss=$(res.fitness)  $param_str")
    end
    println(io)
    println(io, "===== Best result =====")
    println(io, "  Loss = ", best_result.fitness)
    for (i, fp) in enumerate(ModelDef.FITTED_PARAMS)
        println(io, "  $(fp.name) = $(p_best[i])  ($(fp.description))")
    end
    println(io)
    println(io, "===== Information Criteria =====")
    println(io, "  n_obs = $n_obs")
    println(io, "  k     = $k  (fitted params)")
    println(io, "  RSS   = $rss_best")
    println(io, "  logL  = $logL")
    println(io, "  AIC   = $AIC")
    println(io, "  AICc  = $AICc")
    println(io, "  BIC   = $BIC")
end
println("Log saved to $log_path")