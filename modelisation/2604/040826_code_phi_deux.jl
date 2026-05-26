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


## ===== Load data =====
println("\n[ Step 1 ] Load data")
t_step = now()
all_t_H, all_H, all_t_V, all_V = Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}()
dir = joinpath(@__DIR__, "input/xp_input_20°")
n_loaded_fit = 0
for data in ["host","virus"], rep in ["A","B","C"]
    file = "$(data)Data_coevoCondition_Temperature20_Replicate$(rep)_cycle1.csv"
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
all_t_H_rest  = [t[n_change:end] for t in all_t_H]
all_H_rest    = [x[n_change:end] for x in all_H]
all_t_V_first = [t[1:n_change] for t in all_t_V]
all_V_first   = [x[1:n_change] for x in all_V]
all_t_V_rest  = [t[n_change:end] for t in all_t_V]
all_V_rest    = [x[n_change:end] for x in all_V]

all_H0_first = first.(all_H_first)
all_V0_first = first.(all_V_first)
println("  Loaded $n_loaded_fit fitting data files")
log_step("Load fitting data", t_step)


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

function obj_fit1(logθ::AbstractVector)
    try
        p_vec = reconstruct(logθ, model1)
        total_err = 0.0

        for i in eachindex(all_H_first)
            Y0 = log.([all_H0_first[i], 1e-6, 1e-6, all_V0_first[i], 1e-6, 1e-6, 1e-6])
            tspan = (all_t_H_first[i][1], all_t_H_first[i][end])
            prob = ODEProblem(model1.dynamics!, Y0, tspan, p_vec)

            sol = solve(prob, Rodas5(); saveat=all_t_H_first[i])

            H_pred = [logsumexp(sol(t)[1:3]) for t in all_t_H_first[i]]
            V_pred = [logsumexp(sol(t)[4:6]) for t in all_t_H_first[i]]

            total_err += sum((H_pred .- log.(all_H_first[i])).^2)
            total_err += sum((V_pred .- log.(all_V_first[i])).^2)
        end
        return total_err

    catch err
        println("ERROR inside objective: ", err)
        return 1e20
    end
end

res1 = optimize(obj_fit1, lb_fit1, ub_fit1, θ0_fit1, Fminbox(BFGS()))
θ_opt1 = exp.(Optim.minimizer(res1))
println("  ✓ Fit 1 results: μ=$(θ_opt1[1]), φi=$(θ_opt1[2]), β=$(θ_opt1[3]), δ=$(θ_opt1[4]), η=$(θ_opt1[5])")
log_step("Fit 1", t_step)


## ===== Simulate fit 1 =====
println("\n[ Step 3 ] Simulate fit 1")
t_step = now()

t_start1 = all_t_H_first[1][1]
t_end1 = sol_fit1.t[end]
t_fit1 = range(t_start1, t_end1, length=500)

Y0_fit1 = log.([mean(first.(all_H_first)), 1e-6, 1e-6, mean(first.(all_V_first)), 1e-6, 1e-6, 1e-6])
p_fit1 = reconstruct(log.(θ_opt1), model1)

prob_fit1 = ODEProblem(model1.dynamics!, Y0_fit1, (t_start1, t_end1), p_fit1)
sol_fit1 = solve(prob_fit1, Rodas5(); saveat=t_fit1)


## ===== Fit 2 (φi only) =====
println("\n[ Step 4 ] Fit 2 (φi only)")
t_step = now()
fixed_params = Dict(:μ=>θ_opt1[1], :β=>θ_opt1[3], :δ=>θ_opt1[4], :η=>θ_opt1[5])
θ0_fit2 = log.([1.0e-8])
lb_fit2 = log.([1.0e-15])
ub_fit2 = log.([1.0e-5])

Y0_fit2 = sol_fit1.u[end]

function obj_fit2(logθ::AbstractVector)
    try
        p_vec = reconstruct(logθ, model2, fixed_params)
        total_err = 0.0

        for i in eachindex(all_H_rest)
            Y0 = Y0_fit2
            tspan = (all_t_H_rest[i][1], all_t_H_rest[i][end])
            prob = ODEProblem(model2.dynamics!, Y0, tspan, p_vec)

            sol = solve(prob, Rodas5(); saveat=all_t_H_rest[i])

            H_pred = [logsumexp(sol(t)[1:3]) for t in all_t_H_rest[i]]
            V_pred = [logsumexp(sol(t)[4:6]) for t in all_t_H_rest[i]]

            total_err += sum((H_pred .- log.(all_H_rest[i])).^2)
            total_err += sum((V_pred .- log.(all_V_rest[i])).^2)
        end
        return total_err

    catch err
        println("ERROR inside objective: ", err)
        return 1e20
    end
end

res2 = optimize(obj_fit2, lb_fit2, ub_fit2, θ0_fit2, Fminbox(BFGS()))
φi_opt2 = exp.(Optim.minimizer(res2))[1]
println("  ✓ Fit 2 result: φi=$(φi_opt2)")
log_step("Fit 2", t_step)


## ===== Simulate fit 2 =====
println("\n[ Step 5 ] Simulate fit 2")
t_step = now()

t_start2 = all_t_H_rest[1][1]
t_end2 = all_t_H_rest[1][end]
t_fit2 = range(t_start2, t_end2, length=500)

fixed_params_fit2 = Dict(:μ=>θ_opt1[1], :β=>θ_opt1[3], :δ=>θ_opt1[4], :η=>θ_opt1[5])
p_fit2 = reconstruct(log.([φi_opt2]), model2, fixed_params_fit2)

prob_fit2 = ODEProblem(model2.dynamics!, Y0_fit2, (t_start2, t_end2), p_fit2)
sol_fit2 = solve(prob_fit2, Rodas5(); saveat=t_fit2)


## ===== Extraction =====
function extract_all(sol, tvec)
    S = [exp(sol(t)[1]) for t in tvec]
    I = [exp(sol(t)[2]) for t in tvec]
    R = [exp(sol(t)[3]) for t in tvec]
    Vi = [exp(sol(t)[4]) for t in tvec]
    Vdp = [exp(sol(t)[5]) for t in tvec]
    Vdip = [exp(sol(t)[6]) for t in tvec]
    Ev = [exp(sol(t)[7]) for t in tvec]
    return S,I,R,Vi,Vdp,Vdip,Ev
end

S1,I1,R1,Vi1,Vdp1,Vdip1,Ev1 = extract_all(sol_fit1, t_fit1)
S2,I2,R2,Vi2,Vdp2,Vdip2,Ev2 = extract_all(sol_fit2, t_fit2)

# Totaux
y_H_fit1 = S1 .+ I1 .+ R1
y_V_fit1 = Vi1 .+ Vdp1 .+ Vdip1
y_H_fit2 = S2 .+ I2 .+ R2
y_V_fit2 = Vi2 .+ Vdp2 .+ Vdip2

log_step("Simulation", t_step)

## ===== Plot =====
println("\n[ Step 5 ] Plotting")
t_step = now()

host_palette  = cgrad([:darkgreen, :chartreuse])
virus_palette = cgrad([:darkred, :orangered])
col_S = host_palette[0.0]; col_I = host_palette[0.5]; col_R = host_palette[1.0]
col_Vi = virus_palette[0.0]; col_Vdp = virus_palette[0.5]; col_Vdip = virus_palette[1.0]
col_Ev = :blue

p = plot(layout=(1,2), grid=true, yscale=:log10, size=(1700,800),
         xlabel="Time (days)", ylabel="Concentration (parts/mL)",
         legend=:bottomright, margins=15mm)

# ===== Subplot 1 : H & V =====
for i in 1:length(all_H)
    scatter!(p[1], all_t_H[i], all_H[i], color=:green, marker=:circle, alpha=0.5, label=(i==1 ? "H data" : false))
    scatter!(p[1], all_t_V[i], all_V[i], color=:red, marker=:square, alpha=0.5, label=(i==1 ? "V data" : false))
end

plot!(p[1], t_fit1, y_H_fit1, lw=3, color=:green, label="H model cycle 1")
plot!(p[1], t_fit1, y_V_fit1, lw=3, color=:red, label="V model cycle 1")
plot!(p[1], t_fit2, y_H_fit2, lw=3, color=:green, ls=:dash, label="H model cycle 2")
plot!(p[1], t_fit2, y_V_fit2, lw=3, color=:red, ls=:dash, label="V model cycle 2")

title!(p[1], "H and V dynamics")

# ===== Subplot 2 : variables internes =====
plot!(p[2], t_fit1, S1, lw=3, label="S", color=col_S)
plot!(p[2], t_fit1, I1, lw=3, label="I", color=col_I)
plot!(p[2], t_fit1, R1, lw=3, label="R", color=col_R)
plot!(p[2], t_fit1, Vi1, lw=3, label="Vi", color=col_Vi)
plot!(p[2], t_fit1, Vdp1, lw=3, label="Vdp", color=col_Vdp)
plot!(p[2], t_fit1, Vdip1, lw=3, label="Vdip", color=col_Vdip)
plot!(p[2], t_fit1, Ev1, lw=3, label="Ev", color=col_Ev)

plot!(p[2], t_fit2, S2, lw=3, ls=:dash, label=false, color=col_S)
plot!(p[2], t_fit2, I2, lw=3, ls=:dash, label=false, color=col_I)
plot!(p[2], t_fit2, R2, lw=3, ls=:dash, label=false, color=col_R)
plot!(p[2], t_fit2, Vi2, lw=3, ls=:dash, label=false, color=col_Vi)
plot!(p[2], t_fit2, Vdp2, lw=3, ls=:dash, label=false, color=col_Vdp)
plot!(p[2], t_fit2, Vdip2, lw=3, ls=:dash, label=false, color=col_Vdip)
plot!(p[2], t_fit2, Ev2, lw=3, ls=:dash, label=false, color=col_Ev)

title!(p[2], "S, I, R, Vi, Vdp, Vdip, Ev dynamics")

# ===== Y limits globales =====
y_min = 1e-5
all_values = vcat(y_H_fit1, y_V_fit1, y_H_fit2, y_V_fit2,
                 S1,I1,R1,Vi1,Vdp1,Vdip1,Ev1,
                 S2,I2,R2,Vi2,Vdp2,Vdip2,Ev2,
                 vcat(all_H...), vcat(all_V...))

y_min = max(1e-5, minimum(filter(x -> isfinite(x) && x>0, all_values)))
y_max = maximum(filter(x -> isfinite(x) && x>0, all_values))

ylims!(p[1], (y_min, y_max))
ylims!(p[2], (y_min, y_max))

display(p)

run_id = Dates.format(t_global_start, "yyyymmdd-HHMMSS")
plot_path = joinpath(@__DIR__, "080426_output/SIVi_mu/$(run_id)_plot.png")
savefig(p, plot_path)

log_step("Plot", t_step)