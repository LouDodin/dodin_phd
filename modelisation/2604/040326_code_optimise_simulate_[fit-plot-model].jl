using DifferentialEquations
using CSV
using XLSX
using DataFrames
using Plots
using Optim
using Statistics
using Dates
using LogExpFunctions
using Measures
using LaTeXStrings

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
choice_parametrisation = "fixed"                   # "fitted", "fixed" or "fixed_DD"   
choice_fit = 1                                      # 1, 2, 3, 4 or 5 : number of cycles used to fit (if "fitted")
choice_plot = 5                                     # 1, 2, 3, 4 or 5 : number of cycles used for simulation
choice_model = "SIVi"                               # SIVi, SIViVdp, SIRVi_IR

run_id = Dates.format(t_global_start, "yyyymmdd-HHMMSS")
output_dir = joinpath(@__DIR__, "030426_output", "$(choice_parametrisation)_parameters", choice_model)
isdir(output_dir) || mkpath(output_dir)

println("\n" * "="^60)
println("  RUN STARTED: $(Dates.format(t_global_start, "yyyy-mm-dd HH:MM:SS"))")
println("="^60)
println("Choices:")
println("  parametrisation : $choice_parametrisation")
println("  fit cycles      : $choice_fit")
println("  plot cycles     : $choice_plot")
println("  model           : $choice_model")
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


## ===== Load models =====
println("\n[ Step 2 ] Loading model '$choice_model'...")
t_step = now()
include(joinpath(@__DIR__, "models.jl"))
model = MODELS[choice_model]
println("  fit_params  : $(model.fit_params)")
println("  full_params : $(model.full_params)")
log_step("Load model", t_step)


## ===== Extract/reconstruct helpers =====
function extract(dict, params)
    log.([dict[p] for p in params])
end

function reconstruct(logθ, model)
    θ_fit = Dict{Symbol, Float64}()
    for (i,p) in enumerate(model.fit_params)
        θ_fit[p] = exp(logθ[i])
    end
    p_vec = Float64[]
    for p in model.full_params
        if p == :μ;       push!(p_vec, μ)
        elseif p == :k;   push!(p_vec, k)
        else              push!(p_vec, get(θ_fit, p, 0.0))
        end
    end
    return p_vec
end


## ===== Constants and parameters =====
println("\n[ Step 3 ] Setting parameters ($choice_parametrisation)...")
t_step = now()

if choice_parametrisation == "fixed"
    const μ   = 0.4000000000000001;    const k   = 9.784708604680645e7;   const φi  = 1.1156323227682514e-8
    const β   = 69.99999999999996;     const δ   = 0.009999999999999934; const η   = 39.99999999999991
    const εdp = 0.66;   const σdp = 0.26;   const μ_r = 0.7
    const k_r = 1e9;    const ν   = 1;      const α   = 1e-3

elseif choice_parametrisation == "fitted"
    const μ = 0.538001865245161
    const k = 9.784708604680645e7
    const θ_init = Dict(
        :φi => 1.3e-8, :β => 50.0,  :δ => 1e-3,  :η => 10.0,
        :εdp => 1e-2,  :σdp => 1e-3, :μ_r => 0.5, :k_r => 1e8,
        :ν => 1e-12,   :α => 1e-3
    )
    const θ_lower = Dict(
        :φi => 1e-8,  :β => 45.0,  :δ => 1e-4,  :η => 1.0,
        :εdp => 1e-3, :σdp => 1e-6, :μ_r => 1e-3, :k_r => 1e6,
        :ν => 1e-15,  :α => 1e-5
    )
    const θ_upper = Dict(
        :φi => 1e-7,  :β => 70.0,  :δ => 1e-2,  :η => 40.0,
        :εdp => 1e-1, :σdp => 1e-1, :μ_r => 1.0,  :k_r => 1e9,
        :ν => 1e-10,  :α => 1e-1
    )

elseif choice_parametrisation == "fixed_DD"
    const μ   = 0.7;    const k   = 1e9;   const φi  = 2.2e-7
    const β   = 52;     const δ   = 1.5e-2; const η   = 3
    const εdp = 0.66;   const σdp = 0.26;   const μ_r = 0.7
    const k_r = 1e9;    const ν   = 1;      const α   = 1e-3
end

log_step("Set parameters", t_step)


## ===== Model fitting =====
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


if choice_parametrisation == "fitted"

    println("\n[ Step 4 ] Load fitting data ($(choice_fit) cycle(s))...")
    t_step = now()

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
                if isempty(t_data); continue; end
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

    params = model.fit_params
    θ0 = extract(θ_init, params)
    lb = extract(θ_lower, params)
    ub = extract(θ_upper, params)

    println("\n[ Step 5 ] Optimisation...")
    println("  Parameters to fit : $params")
    println("  Initial values    :")
    for (i,p) in enumerate(params)
        println("    $p = $(θ_init[p])  [$(θ_lower[p]), $(θ_upper[p])]")
    end

    println("  Starting optimisation (Fminbox + BFGS)...")
    t_step = now()

    res = optimize(objective_total, lb, ub, θ0, Fminbox(BFGS()))

    t_opt = log_step("Optimisation", t_step)

    p_opt = reconstruct(Optim.minimizer(res), model)
    final_error = Optim.minimum(res)

    println("\n  === Optimisation results ===")
    println("  Converged     : $(Optim.converged(res))")
    println("  Iterations    : $(iteration_count[])")
    println("  Final error   : $(round(final_error, sigdigits=6))")
    println("  Optimized parameters:")
    for (i,p) in enumerate(params)
        v = exp(Optim.minimizer(res)[i])
        println("    $p = $(round(v, sigdigits=5))")
    end
end


## ===== Simulation =====
println("\n[ Step 6 ] Simulation...")
t_step = now()

if choice_parametrisation == "fixed" || choice_parametrisation == "fixed_DD"
    p_sim = Float64[]
    param_map = Dict(:μ=>μ,:k=>k,:φi=>φi,:β=>β,:δ=>δ,:η=>η,
                     :εdp=>εdp,:σdp=>σdp,:μ_r=>μ_r,:k_r=>k_r,:ν=>ν,:α=>α)
    for p in model.full_params
        push!(p_sim, param_map[p])
    end
elseif choice_parametrisation == "fitted"
    p_sim = p_opt
end

t_dilution_all   = [0, 594.5, 834.5, 1051.0, 1363.0, 1603.0] ./ 24
t_dilution   = t_dilution_all[1:choice_plot+1]

H0_mean = [mean([all_H_scatter[3*i+j][1] for j in 1:3]) for i in 0:(length(all_H_scatter)÷3 - 1)]
V0_mean = [mean([all_V_scatter[3*i+j][1] for j in 1:3]) for i in 0:(length(all_V_scatter)÷3 - 1)]

Y0_cycle = log.([H0_mean[1], 1e-6, 1e-6, V0_mean[1], 1e-6, 1e-6, 1e-6])

cycle_props = []

matrix = [[], [], [], [], [], [], [], [], [], []]

for i in 2:length(t_dilution)
    global Y0_cycle
    println("  Simulating cycle $(i-1)/$(length(t_dilution)-1)...")
    t_start = t_dilution[i-1]
    t_end   = t_dilution[i] - 1e-10   # stop just before dilution
    t_cycle = range(t_start, t_end, length=1000)
    prob_cycle = ODEProblem(model.dynamics!, Y0_cycle, (t_start, t_end), p_sim)
    sol_cycle  = solve(prob_cycle, Rodas5(); saveat=t_cycle)

    append!(matrix[1], t_cycle)
    append!(matrix[2], [exp(logsumexp([u[1],u[2],u[3]])) for u in sol_cycle.u])
    append!(matrix[3], exp.(getindex.(sol_cycle.u, 1)))
    append!(matrix[4], exp.(getindex.(sol_cycle.u, 2)))
    append!(matrix[5], exp.(getindex.(sol_cycle.u, 3)))
    append!(matrix[6], [exp(logsumexp([u[4],u[5],u[6]])) for u in sol_cycle.u])
    append!(matrix[7], exp.(getindex.(sol_cycle.u, 4)))
    append!(matrix[8], exp.(getindex.(sol_cycle.u, 5)))
    append!(matrix[9], exp.(getindex.(sol_cycle.u, 6)))
    append!(matrix[10],exp.(getindex.(sol_cycle.u, 7)))

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

        push!(cycle_props, Dict(
            :prop_S=>prop_S,
            :prop_I=>prop_I,
            :prop_R=>prop_R,
            :prop_Vi=>prop_Vi,
            :prop_Vdp=>prop_Vdp,
            :prop_Vdip=>prop_Vdip
        ))

        println("    → Next cycle H0=$(round(H0_next, sigdigits=4)) [S=$(round(S0_next,sigdigits=3)), I=$(round(I0_next,sigdigits=3)), R=$(round(R0_next,sigdigits=3))]")
        println("    → Next cycle V0=$(round(V0_next, sigdigits=4)) [Vi=$(round(Vi0_next,sigdigits=3)), Vdp=$(round(Vdp0_next,sigdigits=3)), Vdip=$(round(Vdip0_next,sigdigits=3))]")
    end
end

log_step("Simulation", t_step)


## ===== Plot =====
println("\n[ Step 7 ] Plot...")
t_step = now()

y_min_raw = 1e-5
all_values = vcat(matrix[2:end]..., all_H_scatter..., all_V_scatter...)
y_max_raw  = maximum(filter(x -> isfinite(x) && x > 0, all_values))
y_lims_raw = (y_min_raw, y_max_raw)
clamp_raw(v) = max.(v, y_min_raw)

exp_min = floor(Int, log10(y_min_raw))
exp_max = ceil(Int,  log10(y_max_raw))
ytick_vals   = [10.0^i for i in exp_min:exp_max]
ytick_labels = [L"10^{%$i}" for i in exp_min:exp_max]

p = plot(
    grid=true, yscale=:log10, yticks=(ytick_vals, ytick_labels),
    layout=(1,2), size=(2000,1000), legend=:bottomright,
    margins=15mm, legendfontsize=14, guidefontsize=14,
    tickfontsize=14, titlefontsize=28,
    plot_title="Model $(choice_model) - Parameters $(choice_parametrisation)",
    plot_titlefontsize=28
)

n_scatter = length(all_H_scatter) ÷ 3
for i in 0:(n_scatter-1)
    for j in 1:3
        scatter!(p[1], all_t_H_scatter[3*i+j], clamp_raw(all_H_scatter[3*i+j]),
                 label=(i==0 && j==1 ? "H flow cytometry" : false), color=:green, marker=:circle, ms=6, alpha=0.6)
        scatter!(p[1], all_t_V_scatter[3*i+j], clamp_raw(all_V_scatter[3*i+j]),
                 label=(i==0 && j==1 ? "V flow cytometry" : false), color=:red, marker=:square, ms=6, alpha=0.6)
    end
end

plot!(p[1], matrix[1], clamp_raw(matrix[2]), label="H model", lw=4, color=:green)
plot!(p[1], matrix[1], clamp_raw(matrix[6]), label="V model", lw=4, color=:red)
title!(p[1], "H and V dynamics")
xlabel!(p[1], "Time (days)")
ylabel!(p[1], "Concentration (parts/mL)")
ylims!(p[1], y_lims_raw)

host_palette  = cgrad([:darkgreen, :chartreuse])
virus_palette = cgrad([:darkred, :orangered])
col_S = host_palette[0.0]; col_I = host_palette[0.5]; col_R = host_palette[1.0]
col_Vi = virus_palette[0.0]; col_Vdp = virus_palette[0.5]; col_Vdip = virus_palette[1.0]
col_Ev = :blue

for (idx, col, lbl) in zip(3:5, [col_S,col_I,col_R], ["S model","I model","R model"])
    plot!(p[2], matrix[1], clamp_raw(matrix[idx]), label=lbl, lw=4, color=col)
end
for (idx, col, lbl) in zip(7:10, [col_Vi,col_Vdp,col_Vdip,col_Ev],
                                  ["Vi model","Vdp model","Vdip model","Ev model"])
    plot!(p[2], matrix[1], clamp_raw(matrix[idx]), label=lbl, lw=4, color=col)
end
title!(p[2], "S, I, R, Vi, Vdp, Vdip and Ev dynamics")
xlabel!(p[2], "Time (days)")
ylabel!(p[2], "Concentration (parts/mL)")
ylims!(p[2], y_lims_raw)

display(p)
log_step("Plot", t_step)


## ===== Save plot and log =====
println("\n[ Step 8 ] Save plot and log...")
t_step = now()

plot_path = joinpath(output_dir, "$(run_id)_fit$(choice_fit)_plot$(choice_plot)_plot.png")
savefig(p, plot_path)
println("  Plot saved to: $plot_path")

t_total = (now() - t_global_start).value / 1000.0

log_path = joinpath(output_dir, "$(run_id)_fit$(choice_fit)_plot$(choice_plot)_log.txt")

open(log_path, "w") do io
    sep = "=" ^ 60

    println(io, sep)
    println(io, "  RUN LOG")
    println(io, "  Date     : $(Dates.format(t_global_start, "yyyy-mm-dd HH:MM:SS"))")
    println(io, "  Run ID   : $(choice_parametrisation)_parameters/$(choice_model)/$(run_id)_fit$(choice_fit)_plot$(choice_plot)")
    println(io, sep)

    println(io, "\n--- CHOICES ---")
    println(io, "  parametrisation : $choice_parametrisation")
    println(io, "  fit cycles      : $choice_fit")
    println(io, "  plot cycles     : $choice_plot")
    println(io, "  model           : $choice_model")

    println(io, "\n--- MODEL INFO ---")
    println(io, "  fit_params  : $(model.fit_params)")
    println(io, "  full_params : $(model.full_params)")

    println(io, "\n--- FIXED PARAMETERS ---")
    println(io, "  μ = $μ")
    println(io, "  k = $k")

    if choice_parametrisation == "fitted"
        println(io, "\n--- INITIAL PARAMETER GUESSES & BOUNDS ---")
        println(io, "  $(rpad("param", 8))  $(rpad("init", 14))  $(rpad("lower", 14))  upper")
        for p in model.fit_params
            println(io, "  $(rpad(string(p),8))  $(rpad(string(θ_init[p]),14))  $(rpad(string(θ_lower[p]),14))  $(θ_upper[p])")
        end

        println(io, "\n--- OPTIMISATION RESULTS ---")
        println(io, "  Converged   : $(Optim.converged(res))")
        println(io, "  Iterations  : $(iteration_count[])")
        println(io, "  Final error : $(Optim.minimum(res))")
        println(io, "\n  Optimized parameter values:")
        println(io, "  $(rpad("param", 8))  $(rpad("fitted", 14))  $(rpad("init", 14))  $(rpad("lower", 14))  upper")
        for (i,p) in enumerate(model.fit_params)
            v = exp(Optim.minimizer(res)[i])
            println(io, "  $(rpad(string(p),8))  $(rpad(string(round(v,sigdigits=6)),14))  $(rpad(string(θ_init[p]),14))  $(rpad(string(θ_lower[p]),14))  $(θ_upper[p])")
        end
        println(io, "\n  Full p_opt vector (order: $(model.full_params)):")
        for (i,p) in enumerate(model.full_params)
            println(io, "    $p = $(p_opt[i])")
        end

    elseif choice_parametrisation in ("fixed", "fixed_DD")
        println(io, "\n--- FIXED PARAMETER VALUES ---")
        for (name, val) in [(:φi,φi),(:β,β),(:δ,δ),(:η,η),(:εdp,εdp),
                             (:σdp,σdp),(:μ_r,μ_r),(:k_r,k_r),(:ν,ν),(:α,α)]
            println(io, "  $name = $val")
        end
    end

    println(io, "\n--- SIMULATION ---")
    println(io, "  Dilution logic  : proportions from end of cycle, applied to experimental H0/V0")
    println(io, "  t_dilution      = $t_dilution")
    println(io, "  Ev fixed at 1e-6 at each cycle start")
    println(io, "")
    println(io, "  $(rpad("cycle", 7))  $(rpad("H0_exp", 12))  $(rpad("prop_S", 10))  $(rpad("prop_I", 10))  $(rpad("prop_R", 10))  $(rpad("V0_exp", 12))  $(rpad("prop_Vi", 10))  $(rpad("prop_Vdp", 10))  prop_Vdip")
    for i in 1:min(choice_plot, length(cycle_props))
        println(io,
            "  $(rpad(string(i), 7))  $(rpad(string(round(H0_mean[i], sigdigits=4)), 12))  $(rpad(string(round(cycle_props[i][:prop_S], sigdigits=4)), 10))  $(rpad(string(round(cycle_props[i][:prop_I], sigdigits=4)), 10))  $(rpad(string(round(cycle_props[i][:prop_R], sigdigits=4)), 10))  $(rpad(string(round(V0_mean[i], sigdigits=4)), 12))  $(rpad(string(round(cycle_props[i][:prop_Vi], sigdigits=4)), 10))  $(rpad(string(round(cycle_props[i][:prop_Vdp], sigdigits=4)), 10))  $(round(cycle_props[i][:prop_Vdip], sigdigits=4))"
        )
    end

    println(io, "\n--- STEP TIMINGS ---")
    for (step, t) in sort(collect(step_times), by=x->x[1])
        println(io, "  $(rpad(step, 30)) : $(round(t, digits=2)) s")
    end
    println(io, "  $(rpad("TOTAL", 30)) : $(round(t_total, digits=2)) s")

    println(io, "\n--- OUTPUT FILES ---")
    println(io, "  plot : $plot_path")
    println(io, "  log  : $log_path")

    println(io, "\n$sep")
    println(io, "  END OF LOG")
    println(io, sep)
end

log_step("Save outputs", t_step)

println("\n" * "="^60)
println("  ALL DONE in $(round(t_total, digits=1))s")
println("  Log saved: $log_path")
println("="^60 * "\n")