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
using DataInterpolations
using Printf
using LaTeXStrings

println(Threads.nthreads())

## ===== Infos =====
color_A = RGB(0.6, 0.8, 1.0)
color_B = RGB(31/255, 119/255, 180/255)
color_C = RGB(0.0, 0.3, 0.7) 
model_color = RGB(255/255, 127/255, 14/255)
data_color = RGB(31/255, 119/255, 180/255)

replicate_colors = [color_A, color_B, color_C]
replicates = ["A", "B", "C"]
n_cycles = 5

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
        df_H = CSV.read(joinpath(@__DIR__, "../input/xp_input_20/hostData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cycle).csv"), DataFrame)
        df_V = CSV.read(joinpath(@__DIR__, "../input/xp_input_20/virusData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cycle).csv"), DataFrame)

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
        u0 = [H[1], V[1], 0]

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

## ===== MODEL =====
function model(dY, Y, p, t)
    ϕ, epsilon = p
    S  = Y[1]
    Vi = Y[2]
    Vd = Y[3]
    dY[1] = r*S*(1 - S/K) - ϕ*S*Vi
    dY[2] = epsilon*β*ϕ*S*Vi - δ*Vi
    dY[3] = (1-epsilon)*β*ϕ*S*Vi - δ*Vd
end

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)

# CHECK THE SIMULATIONS
ϕ = 1e-8
epsilon = 1e-2
p_test  = [ϕ, epsilon]
u0_test = [1E6, 1E7, 0]
prob_test = ODEProblem(model, u0_test, (0.0, 10.0), p_test)
sol_test  = solve(prob_test, Tsit5(), reltol=1e-6, abstol=1e-6, isoutofdomain=isoutofdomain)

pl_sim = plot(layout=(1,2), size=(700,250), margins=5mm)
plot!(pl_sim[1], sol_test.t, sol_test[1,:], label="H", xlabel="time (days)",
      ylabel="abundances (cell/ml)", legend=:bottomright, yscale=:log10, lw=3)
plot!(pl_sim[2], sol_test.t, sol_test[2,:], label="V", xlabel="time (days)",
      ylabel="abundances (virus/ml)", legend=:bottomright, yscale=:log10, lw=3)
#display(pl_sim)


## ===== OBJECTIVE =====
function objective(θ)
    total_err = 0.0
    p_model   = exp.(θ)

    for (key, cyc) in cycles
        (rep, cycle) = key
        if cycle == 1
        # Sélection des points à fitter
        tH = cyc.tH
        H  = cyc.H
        tV = cyc.tV
        V  = cyc.V

        t0     = tH[1]
        t1     = max(tH[end], tV[end])

        prob_cyc = ODEProblem(model, cyc.u0, (t0, t1), p_model)
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
        Vi_pred = [max(sol_cyc(t)[2], 1e-12) for t in tV]
        Vd_pred = [max(sol_cyc(t)[3], 1e-12) for t in tV]

        total_err += sum((log.(S_pred) .- log.(H)).^2) / length(tH)
        total_err += sum((log.(Vi_pred+Vd_pred) .- log.(V)).^2) / length(tV)
        end
    end

    return total_err
end

## ===== MULTI-RUN =====
lower        = log.([1E-15, 1e-10])
upper        = log.([1E-6, 1e-6])
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

pl_fit = plot(
    layout=(2,1),
    size=(1800,1100),

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

# Scatter de toutes les données
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

# Un segment de modèle par cycle, repartant de la moyenne des CI de ce cycle, avec proportion S/R héritée du cycle précédent
global prev_Vi_frac = 1.0  # cycle 1 : 100% infectious
global prev_Vd_frac = 0.0

for cycle in 1:n_cycles
    H0_mean = mean(cycles[(rep, cycle)].H[1] for rep in replicates)
    V0_mean = mean(cycles[(rep, cycle)].V[1] for rep in replicates)
    u0_mean = [H0_mean, V0_mean * prev_Vi_frac, V0_mean * prev_Vd_frac]

    t0 = minimum(cycles[(rep, cycle)].tH[1] for rep in replicates)
    t1 = maximum(max(cycles[(rep, cycle)].tH[end], cycles[(rep, cycle)].tV[end]) for rep in replicates)

    sol_c = solve(
        ODEProblem(model, u0_mean, (t0, t1), exp.(θbest)),
        Rodas5(), reltol=1e-6, abstol=1e-6, isoutofdomain=isoutofdomain
    )

    # Mise à jour des proportions à tf pour le cycle suivant
    Vi_tf = sol_c[2, end]
    Vd_tf = sol_c[3, end]
    total_tf = Vi_tf + Vd_tf
    global prev_Vi_frac = total_tf > 0 ? Vi_tf / total_tf : 1.0
    global prev_Vd_frac = total_tf > 0 ? Vd_tf / total_tf : 0.0

    lbl = cycle == 1 ? "Model" : ""
    plot!(pl_fit[1], sol_c.t, sol_c[1,:], label=lbl, color=model_color, lw=4, alpha=0.7)
    plot!(pl_fit[2], sol_c.t,  sol_c[2,:] .+ sol_c[3,:], label=lbl, color=model_color, lw=4, alpha=0.7)
end

cycle_changes = [24.770833333333332, 34.4375, 43.104166666666664, 55.854166666666664]

for t_change in cycle_changes
    for i in 1:2
        vline!(pl_fit[i], [t_change],
            color=data_color,
            linestyle=:dot,
            lw=2,
            label = t_change == cycle_changes[1] ? "Dilution" : nothing
        )
    end
end

display(pl_fit)
savefig(pl_fit, joinpath(@__DIR__, "020526_output/SViVd_model.png"))

## ===== PLOT phi*Vi/(Vi+Vd) vs phi(t) =====

# Parse du polynôme
n_interv = 3  # à adapter selon ton fichier
poly_file = joinpath(@__DIR__, "240426_output/knots_dilutions_$(n_interv)_polynome_3rep.txt")
function parse_polynome(filepath::String)
    lines = readlines(filepath)

    function extract_float(s)
        m = match(r"([+-])?\s*([0-9]+\.[0-9]+(?:[eE][+-]?[0-9]+)?)", s)
        sign = (m.captures[1] == "-") ? "-" : ""
        return parse(Float64, sign * m.captures[2])
    end

    intervals = []
    i = 1
    while i <= length(lines)
        m = match(r"Interval\s+\[([0-9eE+\-.]+),\s*([0-9eE+\-.]+)\]", lines[i])
        if m !== nothing
            t0 = parse(Float64, m.captures[1])
            a  = extract_float(split(lines[i+1],"=")[2])
            b  = extract_float(split(lines[i+2],"*")[1])
            c  = extract_float(split(lines[i+3],"*")[1])
            d  = extract_float(split(lines[i+4],"*")[1])
            push!(intervals, (t0,a,b,c,d))
            i += 5
        else
            i += 1
        end
    end

    function φ(t)
        for i in 1:length(intervals)-1
            t0 = intervals[i][1]
            t1 = intervals[i+1][1]
            if t ≥ t0 && t < t1
                dt = t - t0
                a,b,c,d = intervals[i][2:end]
                return exp(a + b*dt + c*dt^2 + d*dt^3)
            end
        end
        a,b,c,d = intervals[end][2:end]
        dt = t - intervals[end][1]
        return exp(a + b*dt + c*dt^2 + d*dt^3)
    end

    return φ, intervals
end
φ_global, intervals = parse_polynome(poly_file)

# Couleurs
color_phi_model  = RGB(255/255, 127/255, 14/255)   # orange = modèle
color_phi_import = RGB(0.2, 0.6, 0.2)              # vert = phi importé

pl_phi = plot(
    size=(1800, 550),
    left_margin=15mm,
    right_margin=10mm,
    top_margin=10mm,
    bottom_margin=10mm,
    grid=true,
    xlabel="Time (days)",
    ylabel="Adsorption rate ϕ",
    legend=:topright,
    yscale=:log10,
    xlims=(0, 67),
    ytickfontsize=22,
    legendfontsize=17,
    guidefontsize=22,
    xtickfontsize=22,
)

# Grille de temps globale
t_global = collect(range(0.0, 67.0, length=2000))

# ---- phi importé (polynôme) ----
phi_imported = [φ_global(t) for t in t_global]
plot!(pl_phi, t_global, phi_imported,
    label="ϕ(t) imported",
    color=color_phi_import,
    lw=3,
    linestyle=:dash
)

# ---- phi_eff = phi * Vi/(Vi+Vd) par cycle ----
global prev_Vi_frac = 1.0
global prev_Vd_frac = 0.0

ϕ_fit, _ = exp.(θbest)   # valeur scalaire fittée

for cycle in 1:n_cycles
    H0_mean = mean(cycles[(rep, cycle)].H[1] for rep in replicates)
    V0_mean = mean(cycles[(rep, cycle)].V[1] for rep in replicates)
    u0_mean = [H0_mean, V0_mean * prev_Vi_frac, V0_mean * prev_Vd_frac]

    t0 = minimum(cycles[(rep, cycle)].tH[1] for rep in replicates)
    t1 = maximum(max(cycles[(rep, cycle)].tH[end], cycles[(rep, cycle)].tV[end]) for rep in replicates)

    t_dense = collect(range(t0, t1, length=500))

    sol_c = solve(
        ODEProblem(model, u0_mean, (t0, t1), exp.(θbest)),
        Rodas5(), reltol=1e-6, abstol=1e-6,
        saveat=t_dense,
        isoutofdomain=isoutofdomain
    )

    # Mise à jour fractions pour le cycle suivant
    Vi_tf = sol_c[2, end]
    Vd_tf = sol_c[3, end]
    total_tf = Vi_tf + Vd_tf
    global prev_Vi_frac = total_tf > 0 ? Vi_tf / total_tf : 1.0
    global prev_Vd_frac = total_tf > 0 ? Vd_tf / total_tf : 0.0

    # phi_eff = phi * Vi / (Vi + Vd)
    phi_eff = [
        begin
            Vi = max(sol_c(t)[2], 1e-30)
            Vd = max(sol_c(t)[3], 0.0)
            denom = Vi + Vd
            denom > 0 ? ϕ_fit * Vi / denom : ϕ_fit
        end
        for t in sol_c.t
    ]

    lbl = cycle == 1 ? "ϕ · Vᵢ/(Vᵢ+Vd) model" : ""
    plot!(pl_phi, sol_c.t, phi_eff,
        label=lbl,
        color=color_phi_model,
        lw=3,
        alpha=0.85
    )
end

# Lignes de dilution
for t_change in cycle_changes
    vline!(pl_phi, [t_change],
        color=:gray,
        linestyle=:dot,
        lw=2,
        label = t_change == cycle_changes[1] ? "Dilution" : nothing
    )
end

display(pl_phi)
savefig(pl_phi, joinpath(@__DIR__, "020526_output/phi_eff_vs_phi_imported.png"))