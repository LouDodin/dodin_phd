using CSV
using DataFrames
using DifferentialEquations
using Optim
using Statistics
using Plots
using Measures

# ===== Load data =====
dir_input = joinpath(@__DIR__, "input/xp_input_20°")

t_reps = Vector{Vector{Float64}}()
H_reps = Vector{Vector{Float64}}()

for rep in ("A", "B", "C")
    df = CSV.read(joinpath(dir_input,
        "hostData_hostCondition_Temperature20_Replicate$(rep)_cycle1.csv"),
        DataFrame)

    t = collect(skipmissing(df[:,1])) ./ 24
    h = collect(skipmissing(df[:,2]))

    push!(t_reps, t)
    push!(H_reps, h)
end

# ===== Check alignment =====
@assert all(t_reps[1] == t for t in t_reps) "Temps non alignés"

# ===== Mean trajectory =====
t_data = t_reps[1]
H_data = mean(reduce(hcat, H_reps), dims=2)[:]

# ===== Initial condition =====
S0 = H_data[1]
Y0 = [log(S0)]

# ===== Model =====
function dynamics!(dY, Y, p, t)
    μ, K = p
    S = exp(Y[1])
    dS = μ * S * (1 - S / K)
    dY[1] = dS / S
end

# ===== Loss =====
function loss(logθ)
    μ = exp(logθ[1])
    K = exp(logθ[2])

    prob = ODEProblem(dynamics!, Y0, (minimum(t_data), maximum(t_data)), [μ, K])
    sol = solve(prob, Tsit5(), saveat=t_data)

    pred = exp.(sol[1, :])
    return sum((pred .- H_data).^2)
end

# ===== Fit =====
θ0 = log.([0.5, 1e7])
lb = log.([1e-3, 1e5])
ub = log.([5.0, 1e10])

println("Optimisation en cours...")
res = optimize(loss, lb, ub, θ0, Fminbox(BFGS()))
θ_opt = exp.(Optim.minimizer(res))

println("\n=== FIT RESULTS ===")
println("μ = ", θ_opt[1])
println("K = ", θ_opt[2])

# ===== Simulation =====
prob = ODEProblem(dynamics!, Y0, (minimum(t_data), maximum(t_data)), θ_opt)
sol = solve(prob, Tsit5(), saveat=t_data)
S_pred = exp.(sol[1, :])

# ===== Colors (identiques à ton code) =====
host_palette = cgrad([:darkgreen, :chartreuse])
col_S = host_palette[0.0]

# ===== Global y limits (comme ton code) =====
all_vals = vcat(H_data, S_pred)
y_min = max(1e-5, minimum(filter(x -> isfinite(x) && x > 0, all_vals)))
y_max = maximum(filter(x -> isfinite(x) && x > 0, all_vals))

# ===== Plot =====
fig = plot(layout=(1,2), size=(1600,600),
           xlabel="Time (days)",
           legend=:bottomleft,
           margins=10mm)

# --- Panel 1: H data + fit ---
scatter!(fig[1], t_data, H_data,
    color=:green,
    marker=:circle,
    alpha=0.6,
    label="H data")

plot!(fig[1], t_data, S_pred,
    lw=3,
    color=:green,
    label="Model fit",
    yscale=:log10)

ylabel!(fig[1], "Concentration (parts/mL)")
title!(fig[1], "Host dynamics (mean replicates)")
ylims!(fig[1], (y_min, y_max))

# --- Panel 2: S uniquement (comme ton panel composants) ---
plot!(fig[2], t_data, S_pred,
    lw=3,
    color=col_S,
    label="S",
    yscale=:log10)

scatter!(fig[2], t_data, H_data,
    color=:green,
    marker=:circle,
    alpha=0.4,
    label="data")

ylabel!(fig[2], "Concentration (parts/mL)")
title!(fig[2], "S component")
ylims!(fig[2], (y_min, y_max))

display(fig)