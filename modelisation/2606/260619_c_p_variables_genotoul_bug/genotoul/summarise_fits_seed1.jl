using CSV
using DataFrames
using Printf
using Dates

# ─────────────────────────────────────────────────────────────────────────────
# summarise_fits_c_p.jl
# Parcourt tous les polynomial.txt produits par fit_c_p.jl sous OUTPUT_ROOT
# et génère un tableau récap trié par fitness croissante.
#
# Chaque polynomial.txt encode un fit SR_RS avec :
#   - c(t) et p(t) comme splines cubiques partagées sur n_knots nœuds
#   - ϕ constant (log-espace)
#   - θ = [θ_c (n_knots) ; θ_p (n_knots) ; log(ϕ)]
#
# Usage :
#   julia summarise_fits_c_p.jl                          # utilise le OUTPUT_ROOT par défaut
#   julia summarise_fits_c_p.jl /chemin/vers/output      # chemin custom
# ─────────────────────────────────────────────────────────────────────────────

const OUTPUT_ROOT = length(ARGS) >= 1 ? ARGS[1] :
                    "/work/user/ldodin/modelo/2606/260619/output"

# ── Parseur d'un fichier polynomial.txt (format fit_c_p.jl) ──────────────────

struct RunResult
    combination  :: String       # ex: "0_0_0_0_0"  (issu du nom de dossier nint_...)
    n_interior   :: String       # ex: "[0, 0, 0, 0, 0]"
    total_knots  :: Int
    theta_dim    :: Int          # 2*n_knots + 1
    n_runs       :: Int
    best_fitness :: Float64
    best_seed    :: Int
    export_date  :: String
    # run individuelle
    seed         :: Int
    fitness      :: Float64
    phi          :: Float64      # ϕ = exp(θ[end]) du meilleur candidat
    theta        :: Vector{Float64}
end

function parse_polynomial(path::String)::Vector{RunResult}
    lines = readlines(path)
    results = RunResult[]

    # ── METADATA ──
    export_date  = ""
    n_interior   = ""
    total_knots  = 0
    theta_dim    = 0
    n_runs       = 0
    best_fitness = NaN
    best_seed    = 0

    for l in lines
        if     startswith(l, "export_date")
            export_date  = strip(split(l, "=", limit=2)[end])
        elseif startswith(l, "n_interior/cycle")
            n_interior   = strip(split(l, "=", limit=2)[end])
        elseif startswith(l, "total_knots")
            total_knots  = parse(Int, strip(split(l, "=", limit=2)[end]))
        elseif startswith(l, "θ_dim")
            # ligne :  θ_dim  = 21  (c: 10, p: 10, ϕ: 1)
            m = match(r"θ_dim\s*=\s*(\d+)", l)
            m !== nothing && (theta_dim = parse(Int, m.captures[1]))
        elseif startswith(l, "n_runs")
            n_runs       = parse(Int, strip(split(l, "=", limit=2)[end]))
        elseif startswith(l, "best_fitness")
            best_fitness = parse(Float64, strip(split(l, "=", limit=2)[end]))
        elseif startswith(l, "best_seed")
            best_seed    = parse(Int, strip(split(l, "=", limit=2)[end]))
        end
    end

    # theta_dim peut valoir 0 si la ligne contient θ_dim mais n'a pas matché
    # (caractère Unicode) — on recalcule depuis total_knots
    if theta_dim == 0 && total_knots > 0
        theta_dim = 2 * total_knots + 1
    end

    # Combinaison extraite du chemin (dossier nint_X_X_X_X_X)
    combination = basename(dirname(path))
    combination = replace(combination, "nint_" => "")

    # ── ALL RUNS ──
    # Format : "seed   fitness   [θ_c ; θ_p ; log(ϕ)]"
    in_runs = false
    for l in lines
        if occursin("========== ALL RUNS ==========", l)
            in_runs = true
            continue
        end
        if in_runs && occursin("==========", l)
            in_runs = false
        end
        if !in_runs || isempty(strip(l)) || startswith(strip(l), "seed")
            continue
        end

        m = match(r"^\s*(\d+)\s+([\d.e+\-]+)\s+\[(.+)\]", l)
        m === nothing && continue

        seed    = parse(Int,     m.captures[1])
        fitness = parse(Float64, m.captures[2])
        theta   = parse.(Float64, split(strip(m.captures[3])))

        # ϕ = exp(dernier élément de θ), qui encode log(ϕ)
        phi = isempty(theta) ? NaN : exp(theta[end])

        push!(results, RunResult(
            combination, n_interior, total_knots, theta_dim, n_runs,
            best_fitness, best_seed, export_date,
            seed, fitness, phi, theta
        ))
    end

    # Si aucune run individuelle parsée (vieux format sans section ALL RUNS)
    if isempty(results)
        # Tenter d'extraire ϕ depuis la section PHI
        phi_val = NaN
        for l in lines
            m = match(r"^ϕ\s*=\s*([\d.e+\-]+)", l)
            if m !== nothing
                phi_val = parse(Float64, m.captures[1])
                break
            end
        end
        push!(results, RunResult(
            combination, n_interior, total_knots, theta_dim, n_runs,
            best_fitness, best_seed, export_date,
            best_seed, best_fitness, phi_val, Float64[]
        ))
    end

    return results
end

# ── Collecte de tous les polynomial.txt ──────────────────────────────────────

println("Scanning: $OUTPUT_ROOT")
flush(stdout)

all_files = String[]
for (root, dirs, files) in walkdir(OUTPUT_ROOT)
    for f in files
        if f == "polynomial.txt"
            push!(all_files, joinpath(root, f))
        end
    end
end

println("Found $(length(all_files)) polynomial.txt file(s)")
flush(stdout)

if isempty(all_files)
    println("Aucun fichier trouvé sous $OUTPUT_ROOT")
    exit(1)
end

# ── Parse et assemble ────────────────────────────────────────────────────────

all_results = RunResult[]
for f in all_files
    try
        append!(all_results, parse_polynomial(f))
    catch err
        println("WARNING: failed to parse $f : $err")
        flush(stdout)
    end
end

println("Total runs parsed: $(length(all_results))")
flush(stdout)

# ── Construction du DataFrame ────────────────────────────────────────────────

df = DataFrame(
    rank         = Int[],
    combination  = String[],
    n_interior   = String[],
    total_knots  = Int[],
    theta_dim    = Int[],
    seed         = Int[],
    fitness      = Float64[],
    phi          = Float64[],
    export_date  = String[],
)

for r in sort(all_results, by = x -> x.fitness)
    push!(df, (
        0,
        r.combination,
        r.n_interior,
        r.total_knots,
        r.theta_dim,
        r.seed,
        r.fitness,
        r.phi,
        r.export_date,
    ))
end

df.rank = 1:nrow(df)

# ── Export CSV ───────────────────────────────────────────────────────────────

out_csv = joinpath(OUTPUT_ROOT, "summary_fits_c_p.csv")
CSV.write(out_csv, df)
println("CSV saved → $out_csv")
flush(stdout)

# ── Affichage console ────────────────────────────────────────────────────────

println("\n$(repeat('=', 95))")
println(" SUMMARY — sorted by fitness (best first)  [model: SR_RS, fit_c_p]")
println("$(repeat('=', 95))")
@printf("%-4s  %-15s  %-6s  %-5s  %-5s  %-14s  %-12s\n",
        "rank", "combination", "knots", "θ_dim", "seed", "fitness", "phi")
println(repeat('-', 68))
for row in eachrow(df)
    @printf("%-4d  %-15s  %-6d  %-5d  %-5d  %-14.8e  %-12.4e\n",
            row.rank, row.combination, row.total_knots,
            row.theta_dim, row.seed, row.fitness, row.phi)
end
println(repeat('=', 95))
println("\nDone : $(now())")
flush(stdout)
