using DifferentialEquations
using CSV
using XLSX
using DataFrames
using Plots
using Statistics
using Dates
using LogExpFunctions
using Measures
using LaTeXStrings
using BlackBoxOptim
using Interpolations

gr()


## ===== Timing =====
t_global_start = now()
step_times = Dict{String, Float64}()

function log_step(name::String, t_start)
    elapsed = (now() - t_start).value / 1000.0  # seconds
    step_times[name] = elapsed
    println("  ✓ Done in $(round(elapsed, digits=2))s")
    return elapsed
end


## ===== Choices =====
choice_fit = 1
choice_plot = 5
choice_model = "SIVi_mu"
n_fit = 300

run_id = Dates.format(t_global_start, "yyyymmdd-HHMMSS")
output_dir = joinpath(@__DIR__, "040426_output", choice_model)
isdir(output_dir) || mkpath(output_dir)

println("\n" * "="^60)
println("  RUN STARTED: $(Dates.format(t_global_start, "yyyy-mm-dd HH:MM:SS"))")
println("="^60)
println("Choices:")
println("  number of parametrisation : $n_fit")
println("  fit cycles : $choice_fit")
println("  plot cycles : $choice_plot")
println("  model : $choice_model")
println("="^60 * "\n")


## ===== Load experimental data =====
println("[ Step 1 ] Loading experimental data for plotting...")
t_step = now()

dir = joinpath(@__DIR__, "input/xp_input_20°")
condition = "coevo"
temp = "20"
replicates = ["A","B","C"]

all_t_H_scatter = Vector{Vector{Float64}}()
all_H_scatter   = Vector{Vector{Float64}}()
all_t_V_scatter = Vector{Vector{Float64}}()
all_V_scatter   = Vector{Vector{Float64}}()

n_loaded_scatter = 0
for cycle in 1:choice_plot, rep in replicates
    for data in ["host","virus"]
        file = "$(data)Data_$(condition)Condition_Temperature$(temp)_Replicate$(rep)_Cycle$(cycle).csv"
        path = joinpath(dir, file)
        if isfile(path)
            df = CSV.read(path, DataFrame)
            t_data = collect(skipmissing(df[:,1])) ./ 24
            x_data = collect(skipmissing(df[:,2]))
            if isempty(t_data); continue; end
            if data == "host"
                push!(all_t_H_scatter, t_data)
                push!(all_H_scatter, x_data)
            else
                push!(all_t_V_scatter, t_data)
                push!(all_V_scatter, x_data)
            end
            global n_loaded_scatter += 1
        else
            println("  [WARN] File not found: $file")
        end
    end
end
println("  Loaded $n_loaded_scatter scatter data files ($(length(all_H_scatter)) host, $(length(all_V_scatter)) virus)")
log_step("Load scatter data", t_step)


## ===== Load model =====
println("\n[ Step 2 ] Loading model '$choice_model'...")
t_step = now()
include(joinpath(@__DIR__, "models.jl"))
model = MODELS[choice_model]
println("  fit_params  : $(model.fit_params)")
println("  full_params : $(model.full_params)")
log_step("Load model", t_step)


## ===== Helpers =====
function reconstruct(logθ, model)
    θ_fit = Dict{Symbol, Float64}()
    for (i,p) in enumerate(model.fit_params)
        θ_fit[p] = exp(logθ[i])
    end
    p_vec = Float64[]
    for p in model.full_params
        if p == :k;   push!(p_vec, k)
        else          push!(p_vec, get(θ_fit, p, 0.0))
        end
    end
    return p_vec
end


## ===== Constants and parameters =====
println("\n[ Step 3 ] Parameters...")
t_step = now()

const k   = 9.784708604680645e7

const θ_lower = Dict(:μ=>0.1, :φi=>1e-10, :β=>4.0, :δ=>1e-8, :η=>1.0,
               :εdp=>1e-3, :σdp=>1e-6, :μ_r=>1e-3, :k_r=>1e6,
               :ν=>1e-15, :α=>1e-5)

const θ_upper = Dict(:μ=>0.9, :φi=>1e-6, :β=>100.0, :δ=>1e-2, :η=>100.0,
               :εdp=>1e-1, :σdp=>1e-1, :μ_r=>1.0, :k_r=>1e9,
               :ν=>1e-10, :α=>1e-1)

params = model.fit_params

lower_bounds = [log(θ_lower[p]) for p in params]
upper_bounds = [log(θ_upper[p]) for p in params]

log_step("Parameters", t_step)

## ===== Model fitting =====
println("\n[ Step 4 ] Model fitting on $(choice_fit) cycle(s)...")

# --- Load fitting data ---
println("\n[ Step 4.1 ] Load fitting data...")
t_step = now()

t_max_fit = 186.5

all_t_H = Vector{Vector{Float64}}()
all_H   = Vector{Vector{Float64}}()
all_t_V = Vector{Vector{Float64}}()
all_V   = Vector{Vector{Float64}}()

n_loaded_fit = 0
for cycle in 1:choice_fit, rep in replicates
    for data in ["host","virus"]
        file = "$(data)Data_$(condition)Condition_Temperature$(temp)_Replicate$(rep)_Cycle$(cycle).csv"
        path = joinpath(dir, file)
        if isfile(path)
            df = CSV.read(path, DataFrame)
            t_data = collect(skipmissing(df[:,1])) ./ 24
            x_data = collect(skipmissing(df[:,2]))
            if isempty(t_data)
                continue
            end
            mask = t_data .<= t_max_fit
            t_data = t_data[mask]
            x_data = x_data[mask]
            if isempty(t_data)
                continue
            end
            if data == "host"
                push!(all_t_H, t_data)
                push!(all_H, x_data)
            else
                push!(all_t_V, t_data)
                push!(all_V, x_data)
            end
            global n_loaded_fit += 1
        end
    end
end

all_H0 = first.(all_H)
all_V0 = first.(all_V)

println("  Loaded $n_loaded_fit fitting data files")
log_step("Load fitting data", t_step)

# --- Objective function ---
iteration_count = Ref(0)
last_error      = Ref(Inf)

function objective_total(logθ::AbstractVector)
    try
        p = reconstruct(logθ, model)
        total_err = 0.0

        for i in eachindex(all_H)
            Y0 = log.([all_H0[i], 1e-6, 1e-6, all_V0[i], 1e-6, 1e-6, 1e-6])
            tspan = (all_t_H[i][1], all_t_H[i][end])
            prob  = ODEProblem(model.dynamics!, Y0, tspan, p)

            sol = solve(prob, Rodas5(); saveat=all_t_H[i])

            H_pred = [logsumexp(sol(t)[1:3]) for t in all_t_H[i]]
            V_pred = [logsumexp(sol(t)[4:6]) for t in all_t_H[i]]

            total_err += sum((H_pred .- log.(all_H[i])).^2)
            total_err += sum((V_pred .- log.(all_V[i])).^2)
        end

        last_error[] = total_err

        if iteration_count[] % 50 == 0
            println("    iter $(iteration_count[]) — error = $(round(total_err, sigdigits=5))")
        end

        iteration_count[] += 1
        return total_err

    catch err
        println("ERROR inside objective: ", err)
        return 1e20
    end
end

# --- Optimization ---
println("\n[ Step 4.2 ] Running optimization with Differential Evolution...")
t_step = now()

res_bb = bboptimize(
    objective_total;
    SearchRange = collect(zip(lower_bounds, upper_bounds)),
    NumDimensions = length(params),
    MaxSteps = 2000,
    PopulationSize = n_fit
)

best_logθ = best_candidate(res_bb)
best_err   = best_fitness(res_bb)
p_opt      = reconstruct(best_logθ, model)

println("  Best error = $(round(best_err, digits=4))")
log_step("Optimization", t_step)


## ===== Simulation =====
println("\n[ Step 5 ] Simulation of all fits...")
t_step = now()

t_dilution_all = [0, 594.5, 834.5, 1051.0, 1363.0, 1603.0] ./ 24
t_dilution = t_dilution_all[1:choice_plot+1]

H0_mean = [mean([all_H_scatter[3*i+j][1] for j in 1:3]) for i in 0:(length(all_H_scatter)÷3 - 1)]
V0_mean = [mean([all_V_scatter[3*i+j][1] for j in 1:3]) for i in 0:(length(all_V_scatter)÷3 - 1)]


matrix = []

for (j, logθ) in enumerate(eachcol(res_bb.method_output.population.individuals))

    p_sim = reconstruct(logθ, model)
    Y0_cycle = log.([H0_mean[1], 1e-6, 1e-6, V0_mean[1], 1e-6, 1e-6, 1e-6])

    for i in 2:length(t_dilution)
        println("Simulation : Fit $(j) - Cycle $(i)")
        t_start = t_dilution[i-1]
        t_end   = t_dilution[i] - 1e-10   # stop just before dilution
        t_cycle = range(t_start, t_end, length=1000)
        prob_cycle = ODEProblem(model.dynamics!, Y0_cycle, (t_start, t_end), p_sim)
        sol_cycle = solve(prob_cycle, Rodas5(); saveat=t_cycle)

        push!(matrix, (t_cycle, [exp(logsumexp(u[1:3])) for u in sol_cycle.u], [exp(logsumexp(u[4:6])) for u in sol_cycle.u]))

        # --- Dilution logic ---
        if i < length(t_dilution)
            u_end = sol_cycle.u[end]
            S_end  = exp(u_end[1]); I_end  = exp(u_end[2]); R_end   = exp(u_end[3])
            Vi_end = exp(u_end[4]); Vdp_end= exp(u_end[5]); Vdip_end= exp(u_end[6])

            # Proportions among H
            S_active = S_end >= 2e-6; I_active = I_end >= 2e-6; R_active = R_end >= 2e-6
            if !S_active && !I_active && !R_active
                prop_S, prop_I, prop_R = 1.0, 0.0, 0.0
            else
                H_active = (S_active ? S_end : 0.0) + (I_active ? I_end : 0.0) + (R_active ? R_end : 0.0)
                prop_S = S_active ? S_end / H_active : 0.0
                prop_I = I_active ? I_end / H_active : 0.0
                prop_R = R_active ? R_end / H_active : 0.0
            end

            # Proportions among V
            Vi_active = Vi_end >= 2e-6; Vdp_active = Vdp_end >= 2e-6; Vdip_active = Vdip_end >= 2e-6
            if !Vi_active && !Vdp_active && !Vdip_active
                prop_Vi, prop_Vdp, prop_Vdip = 1.0, 0.0, 0.0
            else
                V_active = (Vi_active ? Vi_end : 0.0) + (Vdp_active ? Vdp_end : 0.0) + (Vdip_active ? Vdip_end : 0.0)
                prop_Vi   = Vi_active   ? Vi_end   / V_active : 0.0
                prop_Vdp  = Vdp_active  ? Vdp_end  / V_active : 0.0
                prop_Vdip = Vdip_active ? Vdip_end / V_active : 0.0
            end

            # New initial conditions from experimental means
            H0_next = H0_mean[i]
            V0_next = V0_mean[i]

            S0_next   = prop_S   > 0.0 ? prop_S   * H0_next : 1e-6
            I0_next   = prop_I   > 0.0 ? prop_I   * H0_next : 1e-6
            R0_next   = prop_R   > 0.0 ? prop_R   * H0_next : 1e-6
            Vi0_next  = prop_Vi  > 0.0 ? prop_Vi  * V0_next : 1e-6
            Vdp0_next = prop_Vdp > 0.0 ? prop_Vdp * V0_next : 1e-6
            Vdip0_next= prop_Vdip> 0.0 ? prop_Vdip* V0_next : 1e-6

            Y0_cycle = log.([S0_next, I0_next, R0_next, Vi0_next, Vdp0_next, Vdip0_next, 1e-6])

            println("    → Next cycle H0=$(round(H0_next, sigdigits=4)) [S=$(round(S0_next,sigdigits=3)), I=$(round(I0_next,sigdigits=3)), R=$(round(R0_next,sigdigits=3))]")
            println("    → Next cycle V0=$(round(V0_next, sigdigits=4)) [Vi=$(round(Vi0_next,sigdigits=3)), Vdp=$(round(Vdp0_next,sigdigits=3)), Vdip=$(round(Vdip0_next,sigdigits=3))]")
        end
    end
end

log_step("Simulation", t_step)


## ===== Combined Plot =====
println("\n[ Step 6 ] Plot...")
t_step = now()

n_sweep = 30
n_rows = 1 + length(params)
n_cols = 2

p = plot(layout=(n_rows, n_cols), size=(1400, 500*n_rows), yscale=:log10, grid=true, margins=10mm)

# --- Row 1: Global fit ---
n_simulated = length(matrix)
palette = cgrad(:viridis, n_simulated)

for j in 1:n_simulated
    t_sim, H_sim, V_sim = matrix[j]
    plot!(p[1,1], t_sim, H_sim, color=:black, lw=1, alpha=0.1, label=false)
    plot!(p[1,2], t_sim, V_sim, color=:black, lw=1, alpha=0.1, label=false)
end
for i in 1:length(all_H_scatter)
    scatter!(p[1,1], all_t_H_scatter[i], all_H_scatter[i], color=:green, label=false)
end
for i in 1:length(all_V_scatter)
    scatter!(p[1,2], all_t_V_scatter[i], all_V_scatter[i], color=:red, label=false)
end
title!(p[1,1], "H — All parameters varying")
title!(p[1,2], "V — All parameters varying")
xlabel!(p, "Time (day)")
ylabel!(p, "Abundance (part/mL)")

# --- Other rows: One parameter at a time ---
θ_opt = Dict(
    :μ  => 0.538001865245161,
    :φi => 1.0887472693771001e-8,
    :β  => 69.99999999999996,
    :δ  => 0.00010000000000000026,
    :η  => 39.99999999999998
)

for (row_idx, param) in enumerate(params)
    p_min = θ_lower[param]
    p_max = θ_upper[param]
    sweep_vals = exp.(range(log(p_min), log(p_max), length=n_sweep))
    colors = cgrad([:yellow, :red], n_sweep)
    row = row_idx + 1
    for val_idx in 1:length(sweep_vals)
        val = sweep_vals[val_idx]
        θ_current = copy(θ_opt)
        θ_current[param] = val
        logθ = [log(θ_current[p]) for p in params]
        p_vec = reconstruct(logθ, model)
        global Y0_cycle = log.([H0_mean[1], 1e-6, 1e-6, V0_mean[1], 1e-6, 1e-6, 1e-6])
        for i in 1:choice_plot
            t_start = t_dilution[i]
            t_end   = t_dilution[i+1] - 1e-10
            t_cycle = range(t_start, t_end, length=400)
            prob_cycle = ODEProblem(model.dynamics!, Y0_cycle, (t_start, t_end), p_vec)
            sol_cycle = try
                solve(prob_cycle, Rodas5(); saveat=t_cycle)
            catch
                continue
            end
            if length(sol_cycle.u) != length(t_cycle) || any(!isfinite, [exp(logsumexp(u[1:3])) for u in sol_cycle.u]) || any(!isfinite, [exp(logsumexp(u[4:6])) for u in sol_cycle.u])
                continue
            end
            H_cycle = [exp(logsumexp(u[1:3])) for u in sol_cycle.u]
            V_cycle = [exp(logsumexp(u[4:6])) for u in sol_cycle.u]
            plot!(p[row,1], t_cycle, H_cycle, color=colors[val_idx], lw=2, alpha=0.8, label=false)
            plot!(p[row,2], t_cycle, V_cycle, color=colors[val_idx], lw=2, alpha=0.8, label=false)

            # dilution
            if i < length(t_dilution)
                u_end = sol_cycle.u[end]
                S_end  = exp(u_end[1]); I_end  = exp(u_end[2]); R_end   = exp(u_end[3])
                Vi_end = exp(u_end[4]); Vdp_end= exp(u_end[5]); Vdip_end= exp(u_end[6])

                # Proportions among H
                S_active = S_end >= 2e-6; I_active = I_end >= 2e-6; R_active = R_end >= 2e-6
                if !S_active && !I_active && !R_active
                    prop_S, prop_I, prop_R = 1.0, 0.0, 0.0
                else
                    H_active = (S_active ? S_end : 0.0) + (I_active ? I_end : 0.0) + (R_active ? R_end : 0.0)
                    prop_S = S_active ? S_end / H_active : 0.0
                    prop_I = I_active ? I_end / H_active : 0.0
                    prop_R = R_active ? R_end / H_active : 0.0
                end

                # Proportions among V
                Vi_active = Vi_end >= 2e-6; Vdp_active = Vdp_end >= 2e-6; Vdip_active = Vdip_end >= 2e-6
                if !Vi_active && !Vdp_active && !Vdip_active
                    prop_Vi, prop_Vdp, prop_Vdip = 1.0, 0.0, 0.0
                else
                    V_active = (Vi_active ? Vi_end : 0.0) + (Vdp_active ? Vdp_end : 0.0) + (Vdip_active ? Vdip_end : 0.0)
                    prop_Vi   = Vi_active   ? Vi_end   / V_active : 0.0
                    prop_Vdp  = Vdp_active  ? Vdp_end  / V_active : 0.0
                    prop_Vdip = Vdip_active ? Vdip_end / V_active : 0.0
                end

                # New initial conditions from experimental means
                H0_next = H0_mean[i]
                V0_next = V0_mean[i]

                S0_next   = prop_S   > 0.0 ? prop_S   * H0_next : 1e-6
                I0_next   = prop_I   > 0.0 ? prop_I   * H0_next : 1e-6
                R0_next   = prop_R   > 0.0 ? prop_R   * H0_next : 1e-6
                Vi0_next  = prop_Vi  > 0.0 ? prop_Vi  * V0_next : 1e-6
                Vdp0_next = prop_Vdp > 0.0 ? prop_Vdp * V0_next : 1e-6
                Vdip0_next= prop_Vdip> 0.0 ? prop_Vdip* V0_next : 1e-6

                Y0_cycle = log.([S0_next, I0_next, R0_next, Vi0_next, Vdp0_next, Vdip0_next, 1e-6])
            end
        end
    end

    # scatter overlay
    for i in 1:length(all_H_scatter)
        scatter!(p[row,1], all_t_H_scatter[i], all_H_scatter[i], color=:green, label=false)
    end
    for i in 1:length(all_V_scatter)
        scatter!(p[row,2], all_t_V_scatter[i], all_V_scatter[i], color=:red, label=false)
    end
    
    title!(p[row,1], "H — $(param) ∈ [$(round(p_min, sigdigits=3)), $(round(p_max, sigdigits=3))]")
    title!(p[row,2], "V — $(param) ∈ [$(round(p_min, sigdigits=3)), $(round(p_max, sigdigits=3))]")
end

display(p)
log_step("Plot", t_step)

## ===== Save plot and log =====
println("\n[ Step 7 ] Save plot and log...")
t_step = now()

plot_path = joinpath(output_dir, "$(run_id)_plot.png")
savefig(p, plot_path)
println("Plot saved to: $plot_path")

t_total = (now() - t_global_start).value / 1000.0

log_path = joinpath(output_dir, "$(run_id)_log.txt")

open(log_path, "w") do io
    println(io, "Run ID: $run_id")
    println(io, "Model: $choice_model")
    println(io, "n_fit: $n_fit")
    println(io, "Fit cycles: $choice_fit")
    println(io, "Plot cycles: $choice_plot")
    println("ALL DONE in $(round(t_total, digits=1))s")
    println(io, "\nStep timings (s):")
    for (k_name, v) in step_times
        println(io, "  $k_name: $(round(v, digits=2))")
    end

    # Paramètres fixes
    println(io, "\nFixed parameters:")
    println(io, "  k = $k")

    # Paramètres fittés avec bornes
    println(io, "\nFitted parameters with bounds:")
    for p in params
        println(io, "  $p: lower=$(θ_lower[p]), upper=$(θ_upper[p])")
    end
end
println("Log saved to: $log_path")