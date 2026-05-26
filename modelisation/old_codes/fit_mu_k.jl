push!(LOAD_PATH, @__DIR__)

using DifferentialEquations
using CSV
using DataFrames
using Plots
using Optim
using Statistics

# ------------------------------
# Dossiers et fichiers
# ------------------------------
dir_input = joinpath(@__DIR__, "C_generated_data")

data = "host"
condition = "host"
temp = "20"
replicates = ["A", "B", "C"]
cycles = ["1", "2", "3", "4", "5", "6", "7"]

# ------------------------------
# Chargement des données (en jours)
# ------------------------------
all_t = Vector{Vector{Float64}}()
all_X = Vector{Vector{Float64}}()
all_rep = String[]
all_cycle = String[]

for rep in replicates
    for cycle in cycles
        input_file = "$(data)Data_$(condition)Condition_Temperature$(temp)_Replicate$(rep)_Cycle$(cycle).csv"
        full_path = joinpath(dir_input, input_file)

        if isfile(full_path)
            df = CSV.read(full_path, DataFrame)

            t_vec = collect(skipmissing(df[:, 1])) ./ 24  # heures → jours
            X_vec = collect(skipmissing(df[:, 2]))

            if !isempty(t_vec) && !isempty(X_vec)
                push!(all_t, t_vec)
                push!(all_X, X_vec)
                push!(all_rep, rep)
                push!(all_cycle, cycle)
            else
                println("Empty data: $full_path")
            end
        else
            println("File not found: $full_path")
        end
    end
end

println("Loaded $(length(all_X)) datasets")

all_X0 = [X[1] for X in all_X]

# ------------------------------
# Modèle de croissance logistique
# ------------------------------
function logistic_growth(X, θ, t)
    mu, k = θ
    mu * X * (1 - X / k)
end

# ------------------------------
# Fonction objectif
# ------------------------------
function objective(θ, model, all_X0, all_t, all_X)
    error = 0.0
    for i in eachindex(all_X)
        prob = ODEProblem(model, all_X0[i], (all_t[i][1], all_t[i][end]), θ)
        sol = solve(prob, Tsit5(), saveat=all_t[i])
        error += sum((sol.u .- all_X[i]).^2)
    end
    return error
end

# ------------------------------
# Optimisation
# ------------------------------
θ_init = [0.7, 1e9] # mu en day^-1, k en cell/mL
θ_bounds = [(0.5, 1e5), (0.9, 1e12)]

lower_bounds = [b[1] for b in θ_bounds]
upper_bounds = [b[2] for b in θ_bounds]

res = optimize(
    θ -> objective(θ, logistic_growth, all_X0, all_t, all_X),
    lower_bounds, upper_bounds, θ_init,
    Fminbox(NelderMead())
)

θ_opt = Optim.minimizer(res)
println("Optimized parameters (mu in day^-1): ", θ_opt)

# ------------------------------
# Simulation : 1 modèle par cycle
# ------------------------------
t_model = Dict{String, Vector{Float64}}()
X_model = Dict{String, Vector{Float64}}()

for cycle in cycles
    idx = findall(all_cycle .== cycle)
    isempty(idx) && continue

    X0_mean = mean(all_X0[idx])
    t = all_t[idx[1]]
    t_dense = range(t[1], t[end], length=200)

    prob = ODEProblem(logistic_growth, X0_mean, (t[1], t[end]), θ_opt)
    sol = solve(prob, Tsit5(), saveat=t_dense)

    t_model[cycle] = t_dense
    X_model[cycle] = sol.u
end

# ------------------------------
# Visualisation
# ------------------------------
rep_colors = Dict("A" => :blue, "B" => :green, "C" => :red)
rep_plotted = Dict("A" => false, "B" => false, "C" => false)

plt = plot(layout = (2, 1), grid = true, size = (800, 900))

# ===== SUBPLOT 1 : données + modèles par cycle =====
for i in eachindex(all_X)
    rep = all_rep[i]
    scatter!(
        plt[1],
        all_t[i],
        all_X[i],
        color = rep_colors[rep],
        markersize = 5,
        label = rep_plotted[rep] ? "" : "Rep $rep"
    )
    rep_plotted[rep] = true
end

for (i, cycle) in enumerate(cycles)
    haskey(t_model, cycle) || continue
    plot!(
        plt[1],
        t_model[cycle],
        X_model[cycle],
        lw = 3,
        color = :black,
        label = i == 1 ? "Model" : ""
    )
end

xlabel!(plt[1], "Time (day)")
ylabel!(plt[1], "Host concentration (cell/mL)")

# ===== SUBPLOT 2 : projection long terme =====
idx_cycle1 = findall(all_cycle .== cycles[1])
X0_long = mean(all_X0[idx_cycle1])

t_long = range(0.0, 60.0, length=500)

prob_long = ODEProblem(
    logistic_growth,
    X0_long,
    (t_long[1], t_long[end]),
    θ_opt
)

sol_long = solve(prob_long, Tsit5(), saveat=t_long)

plot!(
    plt[2],
    t_long,
    sol_long.u,
    lw = 3,
    color = :black,
    label = "Long-term projection"
)

xlabel!(plt[2], "Time (day)")
ylabel!(plt[2], "Host concentration (cell/mL)")

display(plt)