using Dates
using CSV
using DataFrames
using DifferentialEquations
using SciMLBase
using OrdinaryDiffEqRosenbrock
using Statistics
using BlackBoxOptim
using Random
using DataInterpolations
using Printf

# ─────────────────────────────────────────────────────────────────────────────
# Unbuffered stdout : indispensable sur SLURM (pas de terminal interactif)
# ─────────────────────────────────────────────────────────────────────────────
function plog(msg::String)
    println(msg)
    flush(stdout)
end


## ===== CONFIGURATION =====

const DATA_DIR   = "/work/user/ldodin/modelo/input/xp_input_20"
const REPLICATES = ["A", "B", "C"]
const N_CYCLES   = 5
const N_RUNS     = 1

# N_INTERIOR passé en arguments CLI : julia fit_phi.jl 0 1 2 0 1
# Un nœud SLURM = une combinaison = un appel julia avec 5 entiers en ARGS
const N_INTERIOR = length(ARGS) == N_CYCLES ?
    parse.(Int, ARGS) :
    [0, 0, 0, 0, 0]

const OUTPUT_DIR = "/work/user/ldodin/modelo/2606/260604/output/nint_$(join(N_INTERIOR, '_'))"

# Paramètres biologiques (fixés)
const r = 0.574619342477644
const K = 6.675449070379925e7
const β = 144.0
const δ = 0.02

plog("========== CONFIGURATION ==========")
plog("N_INTERIOR = $N_INTERIOR")
plog("OUTPUT_DIR = $OUTPUT_DIR")
plog("N_RUNS     = $N_RUNS")
plog("Threads    = $(Threads.nthreads())")


# ── Vérification : combinaison déjà traitée ? ─────────────────────────────
const POLY_PATH_CHECK = joinpath(OUTPUT_DIR, "polynomial.txt")
if isfile(POLY_PATH_CHECK)
    plog("\n[SKIP] Combinaison déjà traitée : $OUTPUT_DIR")
    plog("[SKIP] Fichier trouvé : $POLY_PATH_CHECK")
    plog("[SKIP] Arrêt immédiat.")
    exit(0)
end


## ===== STRUCTURE =====

struct Cycle
    rep   :: String
    index :: Int
    tH    :: Vector{Float64}
    H     :: Vector{Float64}
    tV    :: Vector{Float64}
    V     :: Vector{Float64}
    u0    :: Vector{Float64}
end


## ===== IMPORT =====

function load_replicate(rep::String)::Vector{Cycle}
    cycles     = Cycle[]
    t_offset   = nothing
    t_end_prev = nothing

    for cyc in 1:N_CYCLES
        df_H = CSV.read(joinpath(DATA_DIR,
            "hostData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cyc).csv"),
            DataFrame)
        df_V = CSV.read(joinpath(DATA_DIR,
            "virusData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cyc).csv"),
            DataFrame)

        tH = df_H[:, 1] ./ 24.0
        H  = df_H[:, 2]
        tV = df_V[:, 1] ./ 24.0
        V  = df_V[:, 2]

        if t_offset === nothing
            t_offset = tH[1]
        end
        tH .-= t_offset
        tV .-= t_offset

        if t_end_prev !== nothing
            gap = tH[1] - t_end_prev
            tH .-= gap
            tV .-= gap
        end
        t_end_prev = tH[end]

        push!(cycles, Cycle(rep, cyc, tH, H, tV, V, [H[1], V[1]]))
    end
    return cycles
end

plog("\n========== LOADING DATA ==========")
all_cycles     = vcat([load_replicate(rep) for rep in REPLICATES]...)
cycles_per_rep = Dict(rep => filter(c -> c.rep == rep, all_cycles) for rep in REPLICATES)
plog("Loaded: $(length(all_cycles)) cycles")

for cyc in all_cycles
    any(cyc.H .<= 0) && plog("WARNING: H <= 0  rep=$(cyc.rep) cycle=$(cyc.index)")
    any(cyc.V .<= 0) && plog("WARNING: V <= 0  rep=$(cyc.rep) cycle=$(cyc.index)")
    any(.!isfinite.(cyc.H)) && plog("WARNING: non-finite H  rep=$(cyc.rep) cycle=$(cyc.index)")
    any(.!isfinite.(cyc.V)) && plog("WARNING: non-finite V  rep=$(cyc.rep) cycle=$(cyc.index)")
end


## ===== CONSTRUCTION DES NŒUDS =====

function build_global_knots()::Vector{Float64}
    knots = Float64[]
    for cyc in 1:N_CYCLES
        cyc_ref = cycles_per_rep["A"][cyc]
        t0 = cyc_ref.tH[1]
        t1 = cyc_ref.tH[end]
        push!(knots, t0)
        if N_INTERIOR[cyc] > 0
            inner = collect(range(t0, t1, length=N_INTERIOR[cyc]+2)[2:end-1])
            append!(knots, inner)
        end
        push!(knots, t1)
    end
    return sort(unique(knots))
end

global_knots = build_global_knots()
n_knots      = length(global_knots)
plog("\n========== KNOTS ==========")
plog("n_knots = $n_knots  →  $(round.(global_knots, digits=2))")

for i in 2:length(global_knots)
    if global_knots[i] <= global_knots[i-1]
        plog("ERROR: bad knot order at i=$i: $(global_knots[i-1]) >= $(global_knots[i])")
    end
end
plog("θ dimension: $n_knots")


## ===== MODÈLE ODE =====

function ode_model!(dY, Y, phi_func, t)
    S, Vi = Y[1], Y[2]
    ϕ = phi_func(t)
    dY[1] = r * S * (1 - S/K) - ϕ * S * Vi
    dY[2] = β * ϕ * S * Vi - δ * Vi
end

function integrate_cycle(cyc::Cycle, phi_func)
    t0 = cyc.tH[1]
    t1 = max(cyc.tH[end], cyc.tV[end])
    sv = sort(unique(vcat(cyc.tH, cyc.tV)))

    prob = ODEProblem(ode_model!, cyc.u0, (t0, t1), phi_func)

    try
        sol = solve(prob, Rodas5(),
                    reltol=1e-6, abstol=1e-6,
                    saveat=sv,
                    isoutofdomain=(u, p, t) -> any(x -> x < 0, u))
        return SciMLBase.successful_retcode(sol) ? sol : nothing
    catch err
        return nothing
    end
end


## ===== FONCTION OBJECTIF =====

function objective(θ::Vector{Float64}, t_knots::Vector{Float64})
    spline   = CubicSpline(θ, t_knots)
    phi_func = t -> exp(spline(clamp(t, t_knots[1], t_knots[end])))
    total_err = 0.0

    for cyc in all_cycles
        sol = integrate_cycle(cyc, phi_func)
        sol === nothing && return 1e12

        S_pred = [max(sol(t)[1], 1e-12) for t in cyc.tH]
        V_pred = [max(sol(t)[2], 1e-12) for t in cyc.tV]

        total_err += sum((log.(S_pred) .- log.(cyc.H)).^2) / length(cyc.H)
        total_err += sum((log.(V_pred) .- log.(cyc.V)).^2) / length(cyc.V)
    end
    return total_err
end


## ===== OPTIMISATION =====
# Chaque run reçoit son propre RNG isolé → thread-safe
# Les runs sont lancés séquentiellement : bboptimize exploite lui-même
# les threads disponibles (via BLAS/LAPACK), pas de conflit de RNG global.
# Pour paralléliser entre combinaisons → 1 nœud SLURM par combinaison
# (voir submit_phi.sh, --array=1-N avec un fichier combinations.txt)

function run_DE(seed::Int, t_knots::Vector{Float64})
    n  = length(t_knots)
    lo = fill(log(1e-12), n)
    hi = fill(log(1e-4),  n)

    Random.seed!(seed)   # seed le RNG global avant chaque run

    res = bboptimize(
        θ -> objective(θ, t_knots);
        SearchRange          = collect(zip(lo, hi)),
        NumDimensions        = n,
        Method               = :xnes,
        PopulationSize       = 1000,
        MaxSteps             = 10000,
        DifferentialWeight   = 0.5,
        CrossoverProbability = 0.9,
        TraceMode            = :silent,
    )
    return (fitness=best_fitness(res), θ=best_candidate(res), seed=seed)
end


## ===== EXPORT POLYNÔME =====

function export_polynomial(path, spline, t_knots, θbest, best_fit, best_seed, all_results)
    T    = spline.t
    A_sp = spline.u
    z    = spline.z

    open(path, "w") do io
        @printf(io, "========== METADATA ==========\n")
        @printf(io, "export_date      = %s\n",   string(now()))
        @printf(io, "n_cycles         = %d\n",   N_CYCLES)
        @printf(io, "replicates       = %s (joint fit)\n", join(REPLICATES, ", "))
        @printf(io, "n_interior/cycle = %s\n",   string(N_INTERIOR))
        @printf(io, "total_knots      = %d\n",   length(t_knots))
        @printf(io, "n_runs           = %d\n",   N_RUNS)
        @printf(io, "best_fitness     = %.8e\n", best_fit)
        @printf(io, "best_seed        = %d\n\n", best_seed)

        @printf(io, "========== ALL RUNS ==========\n")
        @printf(io, "%-6s  %-14s  %s\n", "seed", "fitness", "θ (log-space)")
        for res in sort(all_results, by=r -> r.fitness)
            θ_str = join([@sprintf("%.4e", v) for v in res.θ], "  ")
            @printf(io, "%-6d  %-14.8e  [%s]\n", res.seed, res.fitness, θ_str)
        end
        @printf(io, "\n")

        @printf(io, "========== OPTIMAL KNOTS ==========\n")
        @printf(io, "%-4s  %-14s  %-14s  %s\n", "i", "t_i (days)", "log(φ(t_i))", "φ(t_i) (mL/cell/day)")
        for (i, (t, lv)) in enumerate(zip(t_knots, θbest))
            @printf(io, "%-4d  %-14.4f  %-14.6e  %-14.6e\n", i, t, lv, exp(lv))
        end
        @printf(io, "\n")

        @printf(io, "========== PIECEWISE POLYNOMIAL OF log(φ(t)) ==========\n")
        @printf(io, "φ(t) = exp(P(t)),  P(t) C² cubic spline over %d intervals\n\n", length(T)-1)

        for i in 1:(length(T)-1)
            ti, ti1 = T[i], T[i+1]
            h = ti1 - ti
            a = A_sp[i]
            b = (A_sp[i+1] - A_sp[i])/h - h*(2*z[i] + z[i+1])/6
            c = z[i]/2
            d = (z[i+1] - z[i])/(6*h)

            @printf(io, "Interval [%.4f, %.4f] days:\n", ti, ti1)
            @printf(io, "  log(φ(t)) = %.6e\n",                 a)
            @printf(io, "           + %.6e * (t - %.4f)\n",     b, ti)
            @printf(io, "           + %.6e * (t - %.4f)^2\n",   c, ti)
            @printf(io, "           + %.6e * (t - %.4f)^3\n\n", d, ti)
        end
    end
end


## ===== MAIN =====

mkpath(OUTPUT_DIR)

# Tests de sanité sur la fonction objectif
plog("\n========== OBJECTIVE TESTS ==========")
for phi0 in [1e-12, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4]
    θ_test = fill(log(phi0), n_knots)
    try
        f = objective(θ_test, global_knots)
        plog("phi=$phi0  →  fitness = $f")
    catch err
        plog("phi=$phi0  →  CRASH: $err")
    end
end

# ── Optimisation : N_RUNS séquentiels ─────────────────────────────────────
# Séquentiel = thread-safe garanti.
# Le parallélisme inter-combinaisons se fait au niveau SLURM (--array),
# pas au niveau des threads Julia à l'intérieur d'un nœud.
plog("\n========== OPTIMISATION ($(N_RUNS) runs séquentiels) ==========")

results = Vector{NamedTuple{(:fitness, :θ, :seed), Tuple{Float64, Vector{Float64}, Int}}}(undef, N_RUNS)

for i in 1:N_RUNS
    plog("  run $i / $N_RUNS  (seed=$i)  start=$(now())")
    results[i] = run_DE(i, global_knots)
    plog("  run $i done  fitness=$(round(results[i].fitness, sigdigits=6))")
end

# Meilleur résultat
best_idx = argmin(r.fitness for r in results)
best     = results[best_idx]

plog("\n========== RESULTS ==========")
plog("Best fitness = $(round(best.fitness, sigdigits=6))  (seed=$(best.seed))")
for res in sort(collect(results), by=r -> r.fitness)
    plog("  seed=$(res.seed)  fitness=$(round(res.fitness, sigdigits=6))")
end

# Export
spline    = CubicSpline(best.θ, global_knots)
poly_path = joinpath(OUTPUT_DIR, "polynomial.txt")
export_polynomial(poly_path, spline, global_knots, best.θ,
                  best.fitness, best.seed, collect(results))
plog("\nPolynomial saved → $poly_path")
plog("Done : $(now())")
