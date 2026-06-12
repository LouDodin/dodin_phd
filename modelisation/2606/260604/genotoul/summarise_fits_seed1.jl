using CSV
using DataFrames
using Printf
using Dates

# ─────────────────────────────────────────────────────────────────────────────
# summarise_fits.jl
# Parcourt tous les polynomial.txt sous OUTPUT_ROOT et génère un tableau récap
# trié par fitness croissante.
#
# Usage :
#   julia summarise_fits.jl                          # utilise le OUTPUT_ROOT par défaut
#   julia summarise_fits.jl /chemin/vers/output      # chemin custom
# ─────────────────────────────────────────────────────────────────────────────

const OUTPUT_ROOT = length(ARGS) >= 1 ? ARGS[1] :
                    "/work/user/ldodin/modelo/2606/260604/output"

# ── Parseur d'un fichier polynomial.txt ──────────────────────────────────────

struct RunResult
    combination  :: String      # ex: "0_0_0_0_0"
    n_interior   :: String      # ex: "[0, 0, 0, 0, 0]"
    total_knots  :: Int
    n_runs       :: Int
    best_fitness :: Float64
    best_seed    :: Int
    export_date  :: String
    # toutes les runs individuelles
    seed         :: Int
    fitness      :: Float64
    theta        :: Vector{Float64}
end

function parse_polynomial(path::String)::Vector{RunResult}
    lines = readlines(path)
    results = RunResult[]

    # ── METADATA ──
    get_meta(key) = begin
        m = match(Regex("^$key\\s*=\\s*(.+)"), join(lines, "\n"), 1)
        m === nothing ? "" : strip(m.captures[1])
    end

    # Récupérer chaque champ
    export_date  = ""
    n_interior   = ""
    total_knots  = 0
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
        elseif startswith(l, "n_runs")
            n_runs       = parse(Int, strip(split(l, "=", limit=2)[end]))
        elseif startswith(l, "best_fitness")
            best_fitness = parse(Float64, strip(split(l, "=", limit=2)[end]))
        elseif startswith(l, "best_seed")
            best_seed    = parse(Int, strip(split(l, "=", limit=2)[end]))
        end
    end

    # Combinaison extraite du chemin (dossier nint_X_X_X_X_X)
    combination = basename(dirname(path))
    combination = replace(combination, "nint_" => "")

    # ── ALL RUNS ──
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

        # Format : "1       1.19892406e+01  [-1.8164e+01  ...]"
        m = match(r"^\s*(\d+)\s+([\d.e+\-]+)\s+\[(.+)\]", l)
        m === nothing && continue

        seed    = parse(Int,     m.captures[1])
        fitness = parse(Float64, m.captures[2])
        theta   = parse.(Float64, split(strip(m.captures[3])))

        push!(results, RunResult(
            combination, n_interior, total_knots, n_runs,
            best_fitness, best_seed, export_date,
            seed, fitness, theta
        ))
    end

    # Si aucune run individuelle parsée (vieux format), créer une entrée depuis METADATA
    if isempty(results)
        push!(results, RunResult(
            combination, n_interior, total_knots, n_runs,
            best_fitness, best_seed, export_date,
            best_seed, best_fitness, Float64[]
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
filter!(r -> r.seed == 1, all_results)
println("After filter seed=1 : $(length(all_results)) combinations")
flush(stdout)

# ── Construction du DataFrame ────────────────────────────────────────────────

df = DataFrame(
    rank         = Int[],
    combination  = String[],
    n_interior   = String[],
    total_knots  = Int[],
    fitness      = Float64[],
    export_date  = String[],
)

for r in sort(all_results, by = x -> x.fitness)
    push!(df, (
        0,
        r.combination,
        r.n_interior,
        r.total_knots,
        r.fitness,
        r.export_date,
    ))
end

df.rank = 1:nrow(df)

# ── Export CSV ───────────────────────────────────────────────────────────────

out_csv = joinpath(OUTPUT_ROOT, "summary_fits_seed1.csv")
CSV.write(out_csv, df)
println("CSV saved → $out_csv")
flush(stdout)

# ── Affichage console ────────────────────────────────────────────────────────

println("\n$(repeat('=', 80))")
println(" SUMMARY — sorted by fitness (best first)")
println("$(repeat('=', 80))")
@printf("%-4s  %-15s  %-6s  %-14s\n",
        "rank", "combination", "knots", "fitness")
println(repeat('-', 46))
for row in eachrow(df)
    @printf("%-4d  %-15s  %-6d  %-14.8e\n",
            row.rank, row.combination, row.total_knots, row.fitness)
end
println(repeat('=', 80))
println("\nDone : $(now())")
flush(stdout)
