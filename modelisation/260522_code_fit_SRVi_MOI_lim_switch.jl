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
using DataInterpolations
using Printf
using LaTeXStrings
using Sundials

println(Threads.nthreads())


## ===== Infos =====
n_cycles = 5

color_A = RGB(0.6, 0.8, 1.0)
color_B = RGB(31/255, 119/255, 180/255)
color_C = RGB(0.0, 0.3, 0.7) 
model_color = RGB(255/255, 127/255, 14/255)
data_color = RGB(31/255, 119/255, 180/255)
cycle_colors = [
    RGB(0.95, 0.45, 0.45),  # rose
    RGB(0.95, 0.70, 0.30),  # pêche
    RGB(0.40, 0.78, 0.40),  # vert
    RGB(0.35, 0.60, 0.95),  # bleu
    RGB(0.70, 0.45, 0.95),  # lavande
]

replicate_colors = [color_A, color_B, color_C]
replicates = ["A", "B", "C"]


## ===== Input =====
cycles = Dict{Tuple{String,Int}, NamedTuple}()

t_H_all = Dict{String, Vector{Float64}}()
H_all   = Dict{String, Vector{Float64}}()
t_V_all = Dict{String, Vector{Float64}}()
V_all   = Dict{String, Vector{Float64}}()

for rep in replicates
    t_H_rep = Float64[]
    H_rep   = Float64[]
    t_V_rep = Float64[]
    V_rep   = Float64[]

    t_H_prev_end = nothing

    for cycle in 1:n_cycles
        df_H = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cycle).csv"), DataFrame)
        df_V = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cycle).csv"), DataFrame)

        t_H = df_H[:, 1] ./ 24
        H   = df_H[:, 2]
        t_V = df_V[:, 1] ./ 24
        V   = df_V[:, 2]

        if t_H_prev_end !== nothing
            shift = t_H[1] - t_H_prev_end
            t_H = t_H .- shift
            t_V = t_V .- shift
        end

        t_H_prev_end = t_H[end]

        # Conditions initiales pour ce cycle
        u0 = [H[1], 0, V[1]]

        cycles[(rep, cycle)] = (tH=t_H, H=H, tV=t_V, V=V, u0=u0)

        append!(t_H_rep, t_H)
        append!(H_rep,   H)
        append!(t_V_rep, t_V)
        append!(V_rep,   V)
    end

    t_H_all[rep] = t_H_rep
    H_all[rep]   = H_rep
    t_V_all[rep] = t_V_rep
    V_all[rep]   = V_rep
end


## ===== CHECK THE DATA =====
pl_data = plot(layout=(1,2), size=(900,300), margins=5mm, legend=:topright)

for (i, rep) in enumerate(replicates)
    scatter!(pl_data[1], t_H_all[rep], H_all[rep],
        color=replicate_colors[i], label="Rep $rep",
        xlabel="time (days)", ylabel="abundances (cell/ml)",
        yscale=:log10, ylims=(1e2, 1e8))

    scatter!(pl_data[2], t_V_all[rep], V_all[rep],
        color=replicate_colors[i], label="Rep $rep",
        xlabel="time (days)", ylabel="abundances (part/ml)",
        yscale=:log10, ylims=(1e3, 1e10))
end

#display(pl_data)




## ===== Constants =====
r = 0.5592225270686286
K = 7.29695252684594e7
β = 144
δ = 0.02
ϕ = 8e-9
γ = 0.00023003274074878747
MOI_lim = 106.31350088622402


## ===== Model =====
function model!(dY, Y, p, t)

    eta = p[1]
    S, R, Vi = Y[1], Y[2], Y[3]
    MOI = Vi/(S)

    if MOI < MOI_lim
        dY[1] = r*S*(1 - (S+R)/K) - ϕ*S*Vi + eta*R
        dY[2] = r*R*(1 - (S+R)/K) - eta*R
        dY[3] = β*ϕ*S*Vi - δ*Vi
    else
        dY[1] = r*S*(1 - (S+R)/K) - ϕ*S*Vi - γ*S + eta*R
        dY[2] = γ*S + r*R*(1 - (S+R)/K) - eta*R
        dY[3] = β*ϕ*S*Vi - δ*Vi
    end
end

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)


## ===== OBJECTIVE =====
function objective(θ)
    total_err = 0.0
    p_model   = exp.(θ)

    for (key, cyc) in cycles
        (rep, cycle) = key

        #cycle != 1 && continue

        # Sélection des points à fitter
        tH = cyc.tH
        H  = cyc.H
        tV = cyc.tV
        V  = cyc.V

        t0     = tH[1]
        t1     = max(tH[end], tV[end])

        prob_cyc = ODEProblem(model!, cyc.u0, (t0, t1), p_model)
        sol_cyc  = solve(
            prob_cyc, Rodas5(),
            reltol        = 1e-6,
            abstol        = 1e-6,
            saveat        = sort(unique(vcat(tH, tV))),
            isoutofdomain = isoutofdomain,
        )

        if sol_cyc.retcode != ReturnCode.Success ||
        any(u -> any(x -> !isfinite(x) || x < 0, u), sol_cyc.u)
            return 1e12
        end

        S_pred = [max(sol_cyc(t)[1], 1e-12) for t in tH]
        R_pred = [max(sol_cyc(t)[2], 1e-12) for t in tH]
        Vi_pred = [max(sol_cyc(t)[3], 1e-12) for t in tV]
        #Vd_pred = [max(sol_cyc(t)[4], 1e-12) for t in tV]

        total_err += sum((log.(S_pred+R_pred) .- log.(H)).^2) / length(tH)
        #total_err += sum((log.(Vi_pred+Vd_pred) .- log.(V)).^2) / length(tV)
        total_err += sum((log.(Vi_pred) .- log.(V)).^2) / length(tV)
    end

    return total_err
end


## ===== MULTI-RUN =====
lower        = log.([1e-10])
upper        = log.([1e-3])
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
println("Best params      = ", exp.(θbest))


## ===== PLOT FINAL =====
ytick_vals1   = [10.0^i for i in 2:8]
ytick_labels1 = [L"10^{%$i}" for i in 2:8]

ytick_vals2   = [10.0^i for i in 5:10]
ytick_labels2 = [L"10^{%$i}" for i in 5:10]

ytick_vals3   = [10.0^i for i in -3:5]
ytick_labels3 = [L"10^{%$i}" for i in -3:5]

pl_fit = plot(
    layout=(3,1),
    size=(1800,1600),

    left_margin=15mm,
    right_margin=10mm,
    top_margin=10mm,
    bottom_margin=10mm,

    grid=true,
    yscale=:log10,

    xlims=(0, 67),

    ytickfontsize=26,
    legendfontsize=17,
    guidefontsize=24,
    xtickfontsize=24,
    titlefontsize=24,
    xlabel="Time (days)",
    legend=:bottomright
)

# Scatter de toutes les données (subplots 1 et 2 inchangés)
for (i, rep) in enumerate(replicates)
    for cycle in 1:n_cycles
        cyc = cycles[(rep, cycle)]
        lbl = cycle == 1 ? "Replicate $rep" : ""
        scatter!(pl_fit[1], cyc.tH, cyc.H, label=lbl, color=replicate_colors[i], alpha=0.7,
                 ylabel="Host abundance\n(cell/mL)",
                 ylims=(1e2, 3e8), yticks=(ytick_vals1, ytick_labels1),
                 markershape=:circle, markersize=8)
        scatter!(pl_fit[2], cyc.tV, cyc.V, label=lbl, color=replicate_colors[i], alpha=0.7,
                 ylabel="Virus abundance\n(part/mL)",
                 ylims=(1e5, 1e10), yticks=(ytick_vals2, ytick_labels2),
                 markershape=:circle, markersize=8, legend=(0.13, 0.4))
    end
end

global prev_S_frac = 1.0
global prev_R_frac = 0.0

for cycle in 1:n_cycles
    H0_mean = mean(cycles[(rep, cycle)].H[1] for rep in replicates)
    V0_mean = mean(cycles[(rep, cycle)].V[1] for rep in replicates)
    u0_mean = [H0_mean * prev_S_frac, H0_mean * prev_R_frac, V0_mean]

    t0 = minimum(cycles[(rep, cycle)].tH[1] for rep in replicates)
    t1 = maximum(max(cycles[(rep, cycle)].tH[end], cycles[(rep, cycle)].tV[end]) for rep in replicates)

    sol_c = solve(
        ODEProblem(model!, u0_mean, (t0, t1), exp.(θbest)),
        Rodas5(), reltol=1e-6, abstol=1e-6, isoutofdomain=isoutofdomain
    )

    S_tf = sol_c[1, end]
    R_tf = sol_c[2, end]
    total_H_tf = S_tf + R_tf
    global prev_S_frac = total_H_tf > 0 ? S_tf / total_H_tf : 1.0
    global prev_R_frac = total_H_tf > 0 ? R_tf / total_H_tf : 0.0

    lbl     = cycle == 1 ? "Model H" : ""
    lbl_S   = cycle == 1 ? "Model S" : ""
    lbl_R   = cycle == 1 ? "Model R" : ""
    lbl_V   = cycle == 1 ? "Model Vi" : ""
    lbl_MOI = cycle == 1 ? "Model MOI" : ""

    S_sol  = sol_c[1, :]
    R_sol  = sol_c[2, :]
    Vi_sol = sol_c[3, :]
    H_sol  = S_sol .+ R_sol
    MOI_sol = Vi_sol ./ max.(S_sol, 1e-12)

    # --- Subplot 1 : S, R, S+R ---
    plot!(pl_fit[1], sol_c.t, H_sol,
          label=lbl, color=model_color, lw=4, alpha=0.7)
    plot!(pl_fit[1], sol_c.t, max.(S_sol, 1e-12),
          label=lbl_S, color=RGB(0.2, 0.7, 0.3), lw=2.5, alpha=0.8, linestyle=:dash)
    plot!(pl_fit[1], sol_c.t, max.(R_sol, 1e-12),
          label=lbl_R, color=RGB(0.8, 0.2, 0.2), lw=2.5, alpha=0.8, linestyle=:dashdot)

    # --- Subplot 2 : Virus ---
    plot!(pl_fit[2], sol_c.t, Vi_sol,
          label=lbl_V, color=model_color, lw=4, alpha=0.7)

    # --- Subplot 3 : MOI ---
    plot!(pl_fit[3], sol_c.t, max.(MOI_sol, 1e-12),
          label=lbl_MOI, color=model_color, lw=4, alpha=0.7,
          ylabel="MOI = Vi/S",
          ylims=(1e-3, 1e5), yticks=(ytick_vals3, ytick_labels3),
          yscale=:log10, legend=:bottomleft)
end

# --- Ligne MOI_lim ---
hline!(pl_fit[3], [MOI_lim],
    color=:black,
    linestyle=:dash,
    lw=2,
    label="MOI_lim"
)

# Lignes de dilution sur les 3 subplots
cycle_changes = [24.770833333333332, 34.4375, 43.104166666666664, 55.854166666666664]

for t_change in cycle_changes
    for i in 1:3
        vline!(pl_fit[i], [t_change],
            color=data_color,
            linestyle=:dot,
            lw=2,
            label = t_change == cycle_changes[1] ? "Dilution" : nothing
        )
    end
end
 
display(pl_fit)
savefig(pl_fit, joinpath(@__DIR__, "260522_output/fit_SRVi_switch.png"))