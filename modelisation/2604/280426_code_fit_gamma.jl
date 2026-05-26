## ===== Packages =====
using Dates
using CSV
using DataFrames
using Optim
using DifferentialEquations
using DataInterpolations
using LogExpFunctions
using Statistics
using Plots
using Measures
using BlackBoxOptim
using Base.Threads
using Printf
using Sundials

println("Threads available: ", Threads.nthreads())


## ===== Settings =====
const n_intervals = 3
const n_runs      = 10


## ===== Directories =====
const input_dir  = joinpath(@__DIR__, "input/xp_input_20")
const poly_dir   = joinpath(@__DIR__, "240426_output")
const output_dir = joinpath(@__DIR__, "280426_output")
mkpath(output_dir)   # ensure output directory exists before any savefig call


## ===== Functions =====
struct CycleData
    t_H    :: Vector{Float64}
    H      :: Vector{Float64}
    t_V    :: Vector{Float64}
    V      :: Vector{Float64}
    H0_exp :: Float64
    V0_exp :: Float64
end

function load_cycle(tag::String, cycle::Int)
    hfile = joinpath(input_dir,
        "hostData_coevoCondition_Temperature20_ReplicateA_Cycle$(cycle).csv")
    vfile = joinpath(input_dir,
        "virusData_coevoCondition_Temperature20_ReplicateA_Cycle$(cycle).csv")
    dfH = CSV.read(hfile, DataFrame)
    dfV = CSV.read(vfile, DataFrame)
    tH  = dfH[:, 1] ./ 24.0
    H   = dfH[:, 2]
    tV  = dfV[:, 1] ./ 24.0
    V   = dfV[:, 2]
    return float.(tH), float.(H), float.(tV), float.(V)
end

function shift_cycle!(tH_cur, tV_cur, tH_prev_end)
    gap = tH_cur[1] - tH_prev_end
    tH_cur .-= gap
    tV_cur .-= gap
end

function parse_polynome(filepath::String)
    isfile(filepath) || error("Polynomial file not found: $filepath")
    lines = readlines(filepath)

    metadata = Dict{String,Any}()
    for line in lines
        m = match(r"^(n_knots|n_intervalles|best_fitness|best_seed|n_runs)\s*=\s*(.+)$",
                  strip(line))
        m === nothing && continue
        key = m.captures[1]
        val = strip(m.captures[2])
        metadata[key] = key in ("n_knots", "n_intervalles", "best_seed", "n_runs") ?
                         parse(Int, val) : parse(Float64, val)
    end

    # Robust float extractor: keeps sign, digits, decimal point, exponent
    function extract_float(s::AbstractString)::Float64
        m = match(r"([+-]?\s*[0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)", strip(s))
        m === nothing && error("Cannot extract float from: \"$s\"")
        parse(Float64, replace(m.captures[1], " " => ""))
    end

    intervals = Vector{NTuple{6,Float64}}()
    i = 1
    while i <= length(lines)
        m_iv = match(
            r"Interval\s+\[([0-9eE+\-.]+),\s*([0-9eE+\-.]+)\]\s+days:",
            strip(lines[i]))
        if m_iv !== nothing
            # Need at least 4 more lines for coefficients
            if i + 4 > length(lines)
                @warn "Incomplete interval block at line $i — skipping"
                i += 1; continue
            end
            t0 = parse(Float64, m_iv.captures[1])
            t1 = parse(Float64, m_iv.captures[2])
            # Lines i+1..i+4: "a = ...", "b*...", "c*...", "d*..."
            a  = extract_float(split(lines[i+1], "=")[end])
            b  = extract_float(split(lines[i+2], "*")[1])
            c  = extract_float(split(lines[i+3], "*")[1])
            d  = extract_float(split(lines[i+4], "*")[1])
            push!(intervals, (t0, t1, a, b, c, d))
            i += 5; continue
        end
        i += 1
    end

    isempty(intervals) && error("No intervals found in $filepath")
    println("  Parsed $(length(intervals)) intervals from $(basename(filepath))")
    haskey(metadata, "best_fitness") &&
        println("  fitness = $(metadata["best_fitness"])")

    t_lo = intervals[1][1]
    t_hi = intervals[end][2]

    function phi_raw(t::Real)::Float64
        tc = clamp(Float64(t), t_lo, t_hi)
        # Find the interval containing tc
        idx = length(intervals)
        for k in eachindex(intervals)
            if tc <= intervals[k][2]
                idx = k; break
            end
        end
        t0, _, a, b, c, d = intervals[idx]
        dt = tc - t0
        return exp(a + b*dt + c*dt^2 + d*dt^3)
    end

    return phi_raw, metadata
end

function model!(dY, Y, p, t)
    γ      = p[1]
    S, R, V = Y[1], Y[2], Y[3]
    ϕt     = ϕ_interp(clamp(t, t_grid[1], t_grid[end]))
    N      = S + R
    dY[1]  = r * S * (1 - N / K) - ϕt * N * V - γ * S
    dY[2]  = γ * S + r * R * (1 - N / K)
    dY[3]  = β * ϕt * N * V - δ * V
end

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)

# Returns (sol, idx_H, idx_V) or nothing on failure.
# idx_H / idx_V map each observed time to its position in sol.t.
function solve_cycle(cyc::CycleData, Y0::Vector{Float64}, γ::Float64)
    tspan  = (cyc.t_H[1], cyc.t_H[end])
    # Build save-times: union of host and virus observation times
    t_save = sort(unique(vcat(cyc.t_H, cyc.t_V)))

    prob = ODEProblem(model!, Y0, tspan, [γ])

    sol = try
        solve(prob, CVODE_BDF(linear_solver=:Dense);
              reltol       = 1e-8,
              abstol       = 1e-10,
              saveat       = t_save,
              isoutofdomain = isoutofdomain)
    catch e
        @warn "Solver threw exception: $e"
        return nothing
    end

    if sol.retcode != ReturnCode.Success
        return nothing
    end
    if any(u -> any(x -> !isfinite(x) || x < 0, u), sol.u)
        return nothing
    end

    # Map observation times → indices in sol.t using exact equality on the
    # sorted-unique grid (tolerance 1e-10 to be safe against fp drift)
    tol = 1e-10
    function nearest_idx(t_save::Vector{Float64}, t_obs::Vector{Float64})
        idxs = Vector{Int}(undef, length(t_obs))
        for (j, t) in enumerate(t_obs)
            i = searchsortedfirst(t_save, t)
            # snap to nearest of i-1, i
            if i > length(t_save)
                idxs[j] = length(t_save)
            elseif i == 1
                idxs[j] = 1
            else
                idxs[j] = abs(t_save[i]-t) < abs(t_save[i-1]-t) ? i : i-1
            end
        end
        return idxs
    end

    idx_H = nearest_idx(t_save, cyc.t_H)
    idx_V = nearest_idx(t_save, cyc.t_V)
    (idx_H === nothing || idx_V === nothing) && return nothing

    return sol, idx_H, idx_V
end

# Dilution Reset Between Cycles
function dilution_reset(sol, cyc_next::CycleData)::Vector{Float64}
    u_end  = sol.u[end]
    S_end, R_end = u_end[1], u_end[2]
    N_end  = S_end + R_end
    fS     = N_end > 0 ? S_end / N_end : 0.5
    return [fS * cyc_next.H0_exp,
            (1 - fS) * cyc_next.H0_exp,
            cyc_next.V0_exp]
end

# Multi-Cycle Objective (log-space SSR)
function objective(θ::Vector{Float64})::Float64
    γ   = exp(θ[1])
    ssr = 0.0

    Y0 = [cycles[1].H0_exp, 0.0, cycles[1].V0_exp]

    for (k, cyc) in enumerate(cycles)

        result = solve_cycle(cyc, Y0, γ)
        result === nothing && return 1e12

        sol, _, _ = result

        # --- interpolation directe via sol(t) mais vectorisée ---
        H_pred = Vector{Float64}(undef, length(cyc.t_H))
        V_pred = Vector{Float64}(undef, length(cyc.t_V))

        @inbounds for j in eachindex(cyc.t_H)
            t = clamp(cyc.t_H[j], sol.t[1], sol.t[end])
            u = sol(t)
            H_pred[j] = max(u[1] + u[2], 1e-12)
        end

        @inbounds for j in eachindex(cyc.t_V)
            t = clamp(cyc.t_V[j], sol.t[1], sol.t[end])
            u = sol(t)
            V_pred[j] = max(u[3], 1e-12)
        end

        ssr += sum((log.(H_pred) .- log.(cyc.H)).^2) +
               sum((log.(V_pred) .- log.(cyc.V)).^2)

        # --- reset entre cycles ---
        if k < length(cycles)
            Y0 = dilution_reset(sol, cycles[k+1])
        end
    end

    return ssr
end

function reconstruct_all(γ::Float64)
    sols = Vector{Any}(undef, length(cycles))
    Y0   = [cycles[1].H0_exp, 0.0, cycles[1].V0_exp]
    for (k, cyc) in enumerate(cycles)
        result = solve_cycle(cyc, Y0, γ)
        result === nothing && error("Solver failed on cycle $k during reconstruction")
        sol, _, _ = result
        sols[k]   = sol
        if k < length(cycles)
            Y0 = dilution_reset(sol, cycles[k+1])
        end
    end
    return sols
end


## ===== Input =====
t_HA, HA, t_VA, VA = load_cycle("A", 1)
t_HB, HB, t_VB, VB = load_cycle("B", 2)
t_HC, HC, t_VC, VC = load_cycle("C", 3)
t_HD, HD, t_VD, VD = load_cycle("D", 4)
t_HE, HE, t_VE, VE = load_cycle("E", 5)

# Shift each cycle so it continues from where the previous ended
# (remove the gap between last time of previous cycle and first time of current)
shift_cycle!(t_HB, t_VB, t_HA[end])
shift_cycle!(t_HC, t_VC, t_HB[end])
shift_cycle!(t_HD, t_VD, t_HC[end])
shift_cycle!(t_HE, t_VE, t_HD[end])

# Global time + data vectors (used only for ϕ grid span)
const t_H_all = vcat(t_HA, t_HB, t_HC, t_HD, t_HE)
const H_all   = vcat(HA,   HB,   HC,   HD,   HE)
const t_V_all = vcat(t_VA, t_VB, t_VC, t_VD, t_VE)
const V_all   = vcat(VA,   VB,   VC,   VD,   VE)


## ===== Constants =====
const r = 0.5592225270686286
const K = 7.29695252684594e7
const β = 144.0
const δ = 0.02


## ===== Build ϕ Interpolation =====
const poly_file = joinpath(poly_dir,
    "knots_dilutions_$(n_intervals)_polynome_3rep.txt")

phi_raw, phi_meta = parse_polynome(poly_file)

const N_GRID  = 10_000
const t_grid  = collect(range(t_H_all[1], t_H_all[end]; length=N_GRID))
const phi_grid = phi_raw.(t_grid)

# LinearInterpolation from DataInterpolations: LinearInterpolation(u, t)
const ϕ_interp = LinearInterpolation(phi_grid, t_grid)

let t_check = range(t_H_all[1], t_H_all[end]; length=1_000)
    p = plot(collect(t_check), ϕ_interp.(collect(t_check));
             yscale=:log10, xlabel="Time (days)", ylabel="ϕ(t)",
             title="Pre-interpolated ϕ(t)", legend=false, lw=2)
    savefig(p, joinpath(output_dir, "phi_interp_plot.png"))
    println("  ϕ interpolation plot saved.")
end


## ===== Cycle Data =====

const cycles = [
    CycleData(t_HA, HA, t_VA, VA, HA[1], VA[1]),
    CycleData(t_HB, HB, t_VB, VB, HB[1], VB[1]),
    CycleData(t_HC, HC, t_VC, VC, HC[1], VC[1]),
    CycleData(t_HD, HD, t_VD, VD, HD[1], VD[1]),
    CycleData(t_HE, HE, t_VE, VE, HE[1], VE[1]),
]


## ===== Optimisation =====
const search_range = [(log(1e-5), log(1e-1))]

ResultType = NamedTuple{(:fitness, :θ, :seed), Tuple{Float64, Vector{Float64}, Int}}
results    = Vector{ResultType}(undef, n_runs)

# Note: BlackBoxOptim is not thread-safe when sharing global RNG state.
# Sequential loop is used here; switch to @threads only if you verified thread safety
# with your version of BlackBoxOptim.
for i in 1:n_runs
    seed = 1000 + i
    println("Run $i / $n_runs  (seed $seed)")
    res = bboptimize(
        objective;
        SearchRange          = search_range,
        NumDimensions        = length(search_range),
        Method               = :xnes,
        PopulationSize       = 1000,
        MaxSteps             = 10_000,
        DifferentialWeight   = 0.5,
        CrossoverProbability = 0.9,
        TraceMode            = :silent,
        RandomSeed           = seed,
    )
    results[i] = (fitness = best_fitness(res),
                  θ       = best_candidate(res),
                  seed    = seed)
end

best_idx    = argmin(r.fitness for r in results)
best_result = results[best_idx]
γ_best      = exp(best_result.θ[1])

@printf "\nBest fitness = %.6e\n" best_result.fitness
@printf "Best seed    = %d\n"    best_result.seed
@printf "Best γ       = %.6e\n"  γ_best


## ===== Reconstruct All Cycles with Best γ =====
sols = reconstruct_all(γ_best)


## ===== Plot All Cycles =====
const cycle_labels = ["A", "B", "C", "D", "E"]
pl = plot(layout=(2, 5), size=(1600, 600), margins=6mm, dpi=150)

for (k, (cyc, sol, lbl)) in enumerate(zip(cycles, sols, cycle_labels))
    # Host panel (top row)
    scatter!(pl[1, k], cyc.t_H, cyc.H;
             label="data", yscale=:log10, ms=3,
             xlabel="t (days)", ylabel="cells/ml",
             title="Host $lbl")
    plot!(pl[1, k], sol.t, sol[1, :] .+ sol[2, :];
          label="S+R (model)", lw=2)

    # Virus panel (bottom row)
    scatter!(pl[2, k], cyc.t_V, cyc.V;
             label="data", yscale=:log10, ms=3,
             xlabel="t (days)", ylabel="virus/ml",
             title="Virus $lbl")
    plot!(pl[2, k], sol.t, sol[3, :];
          label="V (model)", lw=2)
end

out_png = joinpath(output_dir, "model_SRV_allcycles.png")
savefig(pl, out_png)
println("\nPlot saved → $out_png")
println("Done.")