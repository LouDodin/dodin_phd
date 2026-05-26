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

println(Threads.nthreads())


## ===== Input =====
df_HA = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle1.csv"), DataFrame)
df_VA = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle1.csv"), DataFrame)
df_HB = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle2.csv"), DataFrame)
df_VB = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle2.csv"), DataFrame)
df_HC = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle3.csv"), DataFrame)
df_VC = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle3.csv"), DataFrame)
df_HD = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle4.csv"), DataFrame)
df_VD = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle4.csv"), DataFrame)
df_HE = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle5.csv"), DataFrame)
df_VE = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle5.csv"), DataFrame)

t_HA = df_HA[1:end,1]./24; HA = df_HA[1:end,2]
t_VA = df_VA[1:end,1]./24; VA = df_VA[1:end,2]

t_HB = df_HB[1:end,1]./24; HB = df_HB[1:end,2]
ratio = (t_HB[1]-t_HA[end])
t_HB = t_HB .- ratio
t_VB = df_VB[1:end,1]./24; VB = df_VB[1:end,2]
t_VB = t_VB .- ratio

t_HC = df_HC[1:end,1]./24; HC = df_HC[1:end,2]
ratio = (t_HC[1]-t_HB[end])
t_HC = t_HC .- ratio
t_VC = df_VC[1:end,1]./24; VC = df_VC[1:end,2]
t_VC = t_VC .- ratio

t_HD = df_HD[1:end,1]./24; HD = df_HD[1:end,2]
ratio = (t_HD[1]-t_HC[end])
t_HD = t_HD .- ratio
t_VD = df_VD[1:end,1]./24; VD = df_VD[1:end,2]
t_VD = t_VD .- ratio

t_HE = df_HE[1:end,1]./24; HE = df_HE[1:end,2]
ratio = (t_HE[1]-t_HD[end])
t_HE = t_HE .- ratio
t_VE = df_VE[1:end,1]./24; VE = df_VE[1:end,2]
t_VE = t_VE .- ratio

# Vecteurs globaux (pour les plots)
t_H = vcat(t_HA, t_HB, t_HC, t_HD, t_HE)
H   = vcat(HA, HB, HC, HD, HE)
t_V = vcat(t_VA, t_VB, t_VC, t_VD, t_VE)
V   = vcat(VA, VB, VC, VD, VE)


## ===== CONSTANTS =====
r = 0.574619342477644
K = 6.675449070379925e7
β = 144
δ = 0.02


## ===== CYCLES : données et CI expérimentales =====
# Chaque cycle : (t_H, H, t_V, V, u0)
# u0 = [H[1], V[1]] du cycle, valeurs expérimentales
cycles = [
    (tH=t_HA, H=HA, tV=t_VA, V=VA, u0=[HA[1], VA[1]]),
    (tH=t_HB, H=HB, tV=t_VB, V=VB, u0=[HB[1], VB[1]]),
    (tH=t_HC, H=HC, tV=t_VC, V=VC, u0=[HC[1], VC[1]]),
    (tH=t_HD, H=HD, tV=t_VD, V=VD, u0=[HD[1], VD[1]]),
    (tH=t_HE, H=HE, tV=t_VE, V=VE, u0=[HE[1], VE[1]]),
]


## ===== NŒUDS DE SPLINE =====
t_start = min(t_H[1],   t_V[1])
t_end   = max(t_H[end], t_V[end])

# Nœuds uniformes — espacés pour une spline lisse
knot_spacing = 8.0  # jours
t_knots = collect(t_start:knot_spacing:t_end)
if t_knots[end] < t_end
    push!(t_knots, t_end)
end
t_knots = [0, 162.5, 666.5, 907, 1195, 1411, t_end*24]./24
n_knots = length(t_knots)
println("Nombre de nœuds de spline : ", n_knots)


## ===== MODÈLE =====
include("model_SV_2.jl")


## ===== FONCTION OBJECTIF =====
function objective_spline(θ)
    log_phi_spline = CubicSpline(θ, t_knots)
    t_lo = t_knots[1]
    t_hi = t_knots[end]
    phi_func(t) = exp(log_phi_spline(clamp(t, t_lo, t_hi)))

    total_err = 0.0

    for cyc in cycles
        t0 = cyc.tH[1]
        t1 = max(cyc.tH[end], cyc.tV[end])
        u0 = cyc.u0

        t_save = sort(unique(vcat(cyc.tH, cyc.tV)))

        prob = ODEProblem(model, u0, (t0, t1), phi_func)
        sol  = solve(
            prob,
            Rodas5(),
            reltol  = 1e-6,
            abstol  = 1e-6,
            saveat  = t_save,
        )

        if sol.retcode != :Success || any(u -> any(x -> !isfinite(x) || x < 0, u), sol.u)
            return 1e12
        end

        S_pred = [max(sol(t)[1], 1e-12) for t in cyc.tH]
        V_pred = [max(sol(t)[2], 1e-12) for t in cyc.tV]

        total_err += sum((log.(S_pred) .- log.(cyc.H)).^2)
        total_err += sum((log.(V_pred) .- log.(cyc.V)).^2)
    end

    return total_err
end


## ===== OPTIMISATION =====
lower        = fill(log(1e-15), n_knots)
upper        = fill(log(1e-6),  n_knots)
search_range = [(lower[i], upper[i]) for i in eachindex(lower)]

function run_DE(seed)
    res = bboptimize(
        objective_spline;
        SearchRange          = search_range,
        NumDimensions        = length(search_range),
        Method               = :xnes,

        PopulationSize       = 1000,
        MaxSteps             = 10000,

        DifferentialWeight   = 0.5,
        CrossoverProbability = 0.9,
        TraceMode            = :silent,
        RandomSeed           = seed
    )
    return (
        fitness = best_fitness(res),
        θ       = best_candidate(res)
    )
end

n_runs  = 100
results = Vector{NamedTuple{(:fitness, :θ), Tuple{Float64, Vector{Float64}}}}(undef, n_runs)

Threads.@threads for i in 1:n_runs
    println("Thread ", threadid(), " démarre run ", i)
    results[i] = run_DE(1000 + i)
end

best_idx    = argmin(r.fitness for r in results)
best_result = results[best_idx]
θbest       = best_result.θ

println("Best error = ", best_result.fitness)
println("Best phi at knots = ", exp.(θbest))


## ===== RECONSTRUCTION DE LA SOLUTION PAR CYCLE =====
log_phi_spline_best = CubicSpline(θbest, t_knots)
phi_func_best(t) = exp(log_phi_spline_best(t))

sols = []
for cyc in cycles
    t0   = cyc.tH[1]
    t1   = max(cyc.tH[end], cyc.tV[end])
    u0   = cyc.u0
    prob = ODEProblem(model, u0, (t0, t1), phi_func_best)
    sol  = solve(prob, Rodas5(), reltol=1e-6, abstol=1e-6)
    push!(sols, sol)
end


## ===== PLOTS =====
t_fine    = range(t_start, t_end, length=1000)
phi_curve = [phi_func_best(t) for t in t_fine]

data_color = RGB(31/255, 119/255, 180/255)
model_color = RGB(255/255, 127/255, 14/255)
phi_color = RGB(255/255, 127/255, 14/255)
cycle_color = RGBA(0.5, 0.5, 0.5, 0.4)

ytick_vals1   = [10.0^i for i in 3:1:9]
ytick_labels1 = [L"10^{%$i}" for i in 3:1:9]
ytick_vals2   = [10.0^i for i in 6:1:10]
ytick_labels2 = [L"10^{%$i}" for i in 3:1:10]
ytick_vals3   = [10.0^i for i in -15:-7]
ytick_labels3 = [L"10^{%$i}" for i in -15:-7]

pl = plot(
    layout=(3,1),
    size=(2000,1500),

    left_margin=15mm,
    right_margin=10mm,
    top_margin=15mm,
    bottom_margin=10mm,

    grid=true,

    yscale=:log10,
    xlims=(t_start, t_end+20),

    ytickfontsize = 26,
    legendfontsize=20,
    guidefontsize=24,
    xtickfontsize=24,
    titlefontsize=24,

    xlabel="Time (days)",
    legend=:bottomright
)

# Délimitations des cycles
cycle_starts = [cyc.tH[1] for cyc in cycles[2:end]]

vline!(pl[1], cycle_starts, lw=1, linestyle=:dash, color=cycle_color, label="Dilution")
vline!(pl[2], cycle_starts, lw=1, linestyle=:dash, color=cycle_color, label="Dilution")
vline!(pl[3], cycle_starts, lw=1, linestyle=:dash, color=cycle_color, label="Dilution")

# Data
scatter!(pl[1], t_H, H,
    label="Host data",
    markershape=:circle,
    markersize=8,
    markerstrokewidth=0,
    color=data_color,
    ylabel="Host abundance\n(cell/mL)",
    ylims=(1e3, 1e9),
    yticks=(ytick_vals1, ytick_labels1)
)

scatter!(pl[2], t_V, V,
    label="Virus data",
    markershape=:square,
    markersize=8,
    markerstrokewidth=0,
    color=data_color,
    ylabel="Virus abundance\n(virion/mL)",
    ylims=(1e6, 1e10),
    yticks=(ytick_vals2, ytick_labels2)
)

# Model
for (i, sol) in enumerate(sols)
    plot!(pl[1], sol.t, sol[1,:], lw=6, color=model_color, label= i==1 ? "Host model" : nothing)
    plot!(pl[2], sol.t, sol[2,:], lw=6, color=model_color, label= i==1 ? "Virus model" : nothing)
end


# phi(t)
scatter!(pl[3], t_knots, exp.(θbest),
    label=" Knots",
    markershape=:circle, markersize=8,
    markerstrokewidth=0,
    color=phi_color
)

plot!(pl[3], t_fine, phi_curve,
    label=" φ(t) spline", lw=4,
    ylabel="φ (ml/(cell.day))",
    color=phi_color,
    ylims=(1e-15, 1e-7),
    yticks=(ytick_vals3, ytick_labels3)
)

println("Saving plot...")
savefig(pl, joinpath(@__DIR__, "230426_output/model_SV_spline_phi.png"))


## ===== EXPORT COEFFICIENTS POLYNOMIAUX =====
T = log_phi_spline_best.t
A = log_phi_spline_best.u
z = log_phi_spline_best.z

println("\n========== POLYNOMIAL DESCRIPTION OF φ(t) ==========")
println("φ(t) piecewise cubic on $(length(T)-1) intervals\n")

for i in 1:(length(T)-1)
    ti  = T[i]
    ti1 = T[i+1]
    h   = ti1 - ti

    a = A[i]
    b = (A[i+1] - A[i])/h - h*(2*z[i] + z[i+1])/6
    c = z[i]/2
    d = (z[i+1] - z[i])/(6*h)

    @printf("Interval [%.4f, %.4f] days:\n", ti, ti1)
    @printf("  log(φ(t)) = %.6e\n", a)
    @printf("       + %.6e * (t - %.4f)\n", b, ti)
    @printf("       + %.6e * (t - %.4f)^2\n", c, ti)
    @printf("       + %.6e * (t - %.4f)^3\n\n", d, ti)
end