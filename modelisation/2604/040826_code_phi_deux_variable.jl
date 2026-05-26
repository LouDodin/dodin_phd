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
    elapsed = (now() - t_start).value / 1000.0
    step_times[name] = elapsed
    println("  ✓ Done in $(round(elapsed, digits=2))s")
    return elapsed
end

println("\n" * "="^60)
println("  RUN STARTED: $(Dates.format(t_global_start, "yyyy-mm-dd HH:MM:SS"))")

## ===== Load model =====
println("\nLoading models")
t_step = now()
include(joinpath(@__DIR__, "models.jl"))

model1 = MODELS["SIVi_mu"]
model2 = MODELS["SIVi_mu_2"]

log_step("Load models", t_step)

## ===== Helpers =====
function reconstruct(logθ, model::ModelSpec, fixed_params=Dict{Symbol,Float64}())
    θ_fit = Dict{Symbol, Float64}()
    idx = 1
    for p in model.fit_params
        if haskey(fixed_params, p)
            θ_fit[p] = fixed_params[p]
        else
            θ_fit[p] = exp(logθ[idx])
            idx += 1
        end
    end
    p_vec = Float64[]
    for p in model.full_params
        if p == :k
            push!(p_vec, k)
        elseif haskey(fixed_params, p)
            push!(p_vec, fixed_params[p])
        else
            push!(p_vec, get(θ_fit, p, 0.0))
        end
    end
    return p_vec
end

function extract_all(sol, tvec)
    S    = [exp(sol(t)[1]) for t in tvec]
    I    = [exp(sol(t)[2]) for t in tvec]
    R    = [exp(sol(t)[3]) for t in tvec]
    Vi   = [exp(sol(t)[4]) for t in tvec]
    Vdp  = [exp(sol(t)[5]) for t in tvec]
    Vdip = [exp(sol(t)[6]) for t in tvec]
    Ev   = [exp(sol(t)[7]) for t in tvec]
    return S, I, R, Vi, Vdp, Vdip, Ev
end

## ===== Load fitting data =====
println("\n[ Step 1 ] Loading fitting data")
t_step = now()
all_t_H, all_H, all_t_V, all_V = Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}()
dir = joinpath(@__DIR__, "input/xp_input_20°")
n_loaded_fit = 0
for data in ["host","virus"], rep in ["A","B","C"]
    file = "$(data)Data_coevoCondition_Temperature20_Replicate$(rep)_Cycle1.csv"
    path = joinpath(dir, file)
    if isfile(path)
        df = CSV.read(path, DataFrame)
        t_data = collect(skipmissing(df[:,1])) ./ 24
        x_data = collect(skipmissing(df[:,2]))
        isempty(t_data) && continue
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

n_change = 5

all_t_H_first = [t[1:n_change] for t in all_t_H]
all_H_first   = [x[1:n_change] for x in all_H]
all_t_H_rest  = [t[n_change+1:end] for t in all_t_H]
all_H_rest    = [x[n_change+1:end] for x in all_H]
all_t_V_first = [t[1:n_change] for t in all_t_V]
all_V_first   = [x[1:n_change] for x in all_V]
all_t_V_rest  = [t[n_change+1:end] for t in all_t_V]
all_V_rest    = [x[n_change+1:end] for x in all_V]

all_H0_first = first.(all_H_first)
all_V0_first = first.(all_V_first)
println("  Loaded $n_loaded_fit fitting data files")
log_step("Load fitting data", t_step)

## ===== Fit 1 (all parameters) =====
println("\n[ Step 2 ] Fit 1 (all parameters)")
t_step = now()
const k = 9.784708604680645e7
params_fit1 = [:μ, :φi, :β, :δ, :η]

const θ_init = Dict(
    :μ=>0.5, :φi => 1.3e-8, :β => 50.0,  :δ => 1e-3,  :η => 10.0,
    :εdp => 1e-2,  :σdp => 1e-3, :μ_r => 0.5, :k_r => 1e8,
    :ν => 1e-12,   :α => 1e-3
)

θ0_fit1 = log.([θ_init[p] for p in params_fit1])
lb_fit1 = log.([0.4, 1.0e-8, 45.0, 1e-4, 1.0])
ub_fit1 = log.([0.7, 1.0e-7, 70.0, 1e-2, 40.0])

obj_fit1 = function(logθ)
    p_vec = reconstruct(logθ, model1)
    err = 0.0
    for i in eachindex(all_H_first)
        Y0 = log.([all_H0_first[i], 1e-6, 1e-6, all_V0_first[i], 1e-6, 1e-6, 1e-6])
        tspan = (all_t_H_first[i][1], all_t_H_first[i][end])
        prob = ODEProblem(model1.dynamics!, Y0, tspan, p_vec)
        sol = try
            solve(prob, Rodas5(); saveat=all_t_H_first[i])
        catch
            return 1e20
        end
        for (t, Hobs) in zip(all_t_H_first[i], all_H_first[i])
            Hpred = sum(exp.(sol(t)[1:3]))
            err += (Hpred - Hobs)^2 / Hobs^2
        end
        for (t, Vobs) in zip(all_t_V_first[i], all_V_first[i])
            Vpred = sum(exp.(sol(t)[4:6]))
            err += (Vpred - Vobs)^2 / Vobs^2
        end
    end
    return err
end

res1 = optimize(obj_fit1, lb_fit1, ub_fit1, θ0_fit1, Fminbox(BFGS()))
θ_opt1 = exp.(Optim.minimizer(res1))
println("  ✓ Fit 1 results: μ=$(θ_opt1[1]), φi=$(θ_opt1[2]), β=$(θ_opt1[3]), δ=$(θ_opt1[4]), η=$(θ_opt1[5])")
log_step("Fit 1", t_step)

## ===== Fit 2 (φi only) =====
println("\n[ Step 3 ] Fit 2 (φi only)")
t_step = now()
fixed_params = Dict(:μ=>θ_opt1[1], :β=>θ_opt1[3], :δ=>θ_opt1[4], :η=>θ_opt1[5])
θ0_fit2 = log.([1.0e-8])
lb_fit2 = log.([1.0e-15])
ub_fit2 = log.([1.0e-5])

obj_fit2 = function(logθ)
    p_vec = reconstruct(logθ, model2, fixed_params)
    err = 0.0
    for i in eachindex(all_H_rest)
        Y0 = log.([all_H_rest[i][1], 1e-6, 1e-6, all_V_rest[i][1], 1e-6, 1e-6, 1e-6])
        tspan = (all_t_H_rest[i][1], all_t_H_rest[i][end])
        prob = ODEProblem(model2.dynamics!, Y0, tspan, p_vec)
        sol = try
            solve(prob, Rodas5(); saveat=all_t_H_rest[i])
        catch
            return 1e20
        end
        for (t, Hobs) in zip(all_t_H_rest[i], all_H_rest[i])
            Hpred = sum(exp.(sol(t)[1:3]))
            err += (Hpred - Hobs)^2 / Hobs^2
        end
        for (t, Vobs) in zip(all_t_V_rest[i], all_V_rest[i])
            Vpred = sum(exp.(sol(t)[4:6]))
            err += (Vpred - Vobs)^2 / Vobs^2
        end
    end
    return err
end

res2 = optimize(obj_fit2, lb_fit2, ub_fit2, θ0_fit2, Fminbox(BFGS()))
φi_opt2 = exp.(Optim.minimizer(res2))[1]
println("  ✓ Fit 2 result: φi=$(φi_opt2)")
log_step("Fit 2", t_step)

## ===== Simulation =====
println("\n[ Step 4 ] Simulation")
t_step = now()

# --- Cycle 1 (unchanged, single run with fit 1 parameters) ---
t_start1  = 0.0
t_end1    = 3.7708333333333335 + (6.770833333333333 - 3.7708333333333335) * 0
t_cycle1  = range(t_start1, stop=t_end1, length=500)
Y0_cycle1 = log.([mean(first.(all_H_first)), 1e-6, 1e-6, mean(first.(all_V_first)), 1e-6, 1e-6, 1e-6])
p_cycle1  = reconstruct(log.(θ_opt1), model1)
prob_cycle1 = ODEProblem(model1.dynamics!, Y0_cycle1, (t_start1, t_end1), p_cycle1)
sol_cycle1  = solve(prob_cycle1, Rodas5())

S1, I1, R1, Vi1, Vdp1, Vdip1, Ev1 = extract_all(sol_cycle1, t_cycle1)
y_H_cycle1 = S1 .+ I1 .+ R1
y_V_cycle1 = Vi1 .+ Vdp1 .+ Vdip1

# --- Cycle 2 : φi sweep (50 values, DE bounds lb_fit2..ub_fit2) ---
n_sweep      = 50
φi_sweep     = exp.(range(log(exp(lb_fit2[1])), log(exp(ub_fit2[1])), length=n_sweep))
sweep_colors = cgrad([:yellow, :red], n_sweep)

t_start2 = 6.770833333333333 - (6.770833333333333 - 3.7708333333333335) * 1
t_end2   = 24.770833333333332
t_cycle2 = range(t_start2, stop=t_end2, length=500)

Y0_cycle2_base = sol_cycle1(t_end1)   # handoff point from cycle 1 (log-space)

fixed_params_cycle2 = Dict(:μ=>θ_opt1[1], :β=>θ_opt1[3], :δ=>θ_opt1[4], :η=>θ_opt1[5])

# Storage for sweep trajectories
sweep_H    = Vector{Vector{Float64}}()
sweep_V    = Vector{Vector{Float64}}()
sweep_S    = Vector{Vector{Float64}}()
sweep_I    = Vector{Vector{Float64}}()
sweep_R    = Vector{Vector{Float64}}()
sweep_Vi   = Vector{Vector{Float64}}()
sweep_Vdp  = Vector{Vector{Float64}}()
sweep_Vdip = Vector{Vector{Float64}}()
sweep_Ev   = Vector{Vector{Float64}}()

println("  Running φi sweep ($n_sweep values) for cycle 2...")
for (idx, φi_val) in enumerate(φi_sweep)
    p_c2   = reconstruct(log.([φi_val]), model2, fixed_params_cycle2)
    prob_c2 = ODEProblem(model2.dynamics!, Y0_cycle2_base, (t_start2, t_end2), p_c2)
    sol_c2 = try
        solve(prob_c2, Rodas5(); saveat=t_cycle2)
    catch
        push!(sweep_H,    fill(NaN, length(t_cycle2)))
        push!(sweep_V,    fill(NaN, length(t_cycle2)))
        push!(sweep_S,    fill(NaN, length(t_cycle2)))
        push!(sweep_I,    fill(NaN, length(t_cycle2)))
        push!(sweep_R,    fill(NaN, length(t_cycle2)))
        push!(sweep_Vi,   fill(NaN, length(t_cycle2)))
        push!(sweep_Vdp,  fill(NaN, length(t_cycle2)))
        push!(sweep_Vdip, fill(NaN, length(t_cycle2)))
        push!(sweep_Ev,   fill(NaN, length(t_cycle2)))
        continue
    end
    S2, I2, R2, Vi2, Vdp2, Vdip2, Ev2 = extract_all(sol_c2, t_cycle2)
    push!(sweep_H,    S2 .+ I2 .+ R2)
    push!(sweep_V,    Vi2 .+ Vdp2 .+ Vdip2)
    push!(sweep_S,    S2);    push!(sweep_I,    I2);    push!(sweep_R,    R2)
    push!(sweep_Vi,   Vi2);   push!(sweep_Vdp,  Vdp2);  push!(sweep_Vdip, Vdip2)
    push!(sweep_Ev,   Ev2)
end

log_step("Simulation", t_step)

## ===== Plot =====
println("\n[ Step 5 ] Plotting")
t_step = now()

host_palette  = cgrad([:darkgreen, :chartreuse])
virus_palette = cgrad([:darkred, :orangered])
col_S    = host_palette[0.0];   col_I    = host_palette[0.5];   col_R    = host_palette[1.0]
col_Vi   = virus_palette[0.0];  col_Vdp  = virus_palette[0.5];  col_Vdip = virus_palette[1.0]
col_Ev   = :blue

p = plot(layout=(1,2), grid=true, yscale=:log10, size=(1700, 800),
         xlabel="Time (days)", ylabel="Concentration (parts/mL)",
         legend=:bottomright, margins=15mm)

# ===== Subplot 1 : H & V =====
# Scatter data
for i in 1:length(all_H)
    scatter!(p[1], all_t_H[i], all_H[i], color=:green, marker=:circle, alpha=0.5,
             label=(i==1 ? "H data" : false))
    scatter!(p[1], all_t_V[i], all_V[i], color=:red,   marker=:square, alpha=0.5,
             label=(i==1 ? "V data" : false))
end
# Cycle 1 (single)
plot!(p[1], t_cycle1, y_H_cycle1, lw=3, color=:green, label="H model cycle 1")
plot!(p[1], t_cycle1, y_V_cycle1, lw=3, color=:red,   label="V model cycle 1")
# Cycle 2 sweep (yellow→red)
for idx in 1:n_sweep
    plot!(p[1], collect(t_cycle2), sweep_H[idx], lw=1.5, color=sweep_colors[idx], alpha=0.75, label=false)
    plot!(p[1], collect(t_cycle2), sweep_V[idx], lw=1.5, color=sweep_colors[idx], alpha=0.75, label=false)
end
title!(p[1], "H and V — cycle 2 φi sweep\n(yellow → red = $(round(φi_sweep[1], sigdigits=2)) → $(round(φi_sweep[end], sigdigits=2)))")

# ===== Subplot 2 : internal variables =====
# Cycle 1 (single, solid)
plot!(p[2], t_cycle1, S1,    lw=3, label="S",    color=col_S)
plot!(p[2], t_cycle1, I1,    lw=3, label="I",    color=col_I)
plot!(p[2], t_cycle1, R1,    lw=3, label="R",    color=col_R)
plot!(p[2], t_cycle1, Vi1,   lw=3, label="Vi",   color=col_Vi)
plot!(p[2], t_cycle1, Vdp1,  lw=3, label="Vdp",  color=col_Vdp)
plot!(p[2], t_cycle1, Vdip1, lw=3, label="Vdip", color=col_Vdip)
plot!(p[2], t_cycle1, Ev1,   lw=3, label="Ev",   color=col_Ev)
# Cycle 2 sweep: one colour per φi value, all variables
for idx in 1:n_sweep
    c = sweep_colors[idx]
    plot!(p[2], collect(t_cycle2), sweep_S[idx],    lw=1, color=c, alpha=0.6, label=false)
    plot!(p[2], collect(t_cycle2), sweep_I[idx],    lw=1, color=c, alpha=0.6, label=false)
    plot!(p[2], collect(t_cycle2), sweep_R[idx],    lw=1, color=c, alpha=0.6, label=false)
    plot!(p[2], collect(t_cycle2), sweep_Vi[idx],   lw=1, color=c, alpha=0.6, label=false)
    plot!(p[2], collect(t_cycle2), sweep_Vdp[idx],  lw=1, color=c, alpha=0.6, label=false)
    plot!(p[2], collect(t_cycle2), sweep_Vdip[idx], lw=1, color=c, alpha=0.6, label=false)
    plot!(p[2], collect(t_cycle2), sweep_Ev[idx],   lw=1, color=c, alpha=0.6, label=false)
end
title!(p[2], "S, I, R, Vi, Vdp, Vdip, Ev — cycle 2 φi sweep\n(yellow → red = low → high φi)")

# ===== Global Y limits =====
all_values = vcat(
    y_H_cycle1, y_V_cycle1,
    S1, I1, R1, Vi1, Vdp1, Vdip1, Ev1,
    vcat(sweep_H...), vcat(sweep_V...),
    vcat(sweep_S...), vcat(sweep_I...), vcat(sweep_R...),
    vcat(sweep_Vi...), vcat(sweep_Vdp...), vcat(sweep_Vdip...), vcat(sweep_Ev...),
    vcat(all_H...), vcat(all_V...)
)

y_min = max(1e-5, minimum(filter(x -> isfinite(x) && x > 0, all_values)))
y_max = maximum(filter(x -> isfinite(x) && x > 0, all_values))
ylims!(p[1], (y_min, y_max))
ylims!(p[2], (y_min, y_max))

display(p)

run_id    = Dates.format(t_global_start, "yyyymmdd-HHMMSS")
plot_path = joinpath(@__DIR__, "080426_output/SIVi_mu/$(run_id)_plot.png")
mkpath(dirname(plot_path))
savefig(p, plot_path)
println("  Plot saved: $plot_path")

log_step("Plot", t_step)