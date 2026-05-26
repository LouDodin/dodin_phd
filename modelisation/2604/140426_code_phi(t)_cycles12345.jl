## ===== Package import =====
using Dates
using CSV
using DataFrames
using Optim
using DifferentialEquations
using LogExpFunctions
using Measures
using Plots
using Statistics

t_global_start = now()

const RUN_FIT_1       = false
const RUN_FIT_OTHERS  = true
const FIT_PHI_NEXT_CYCLES = true   # true  → optimise phi à chaque segment
                                   # false → optimise phi cycle 1 seulement, puis garde le dernier phi constant

cycles = [1,2]

## ===== Imports =====

dir_input = joinpath(@__DIR__, "input/xp_input_20°")

t_H = Vector{Vector{Float64}}()
H   = Vector{Vector{Float64}}()
t_V = Vector{Vector{Float64}}()
V   = Vector{Vector{Float64}}()

# dim: [cycle][rep] — on stocke par (cycle, rep) dans l'ordre cycle1A, cycle1B, cycle1C, cycle2A, ...
for cycle in cycles, rep in ("A", "B", "C")
    for data in ("host", "virus")
        df = CSV.read(joinpath(dir_input,
            "$(data)Data_coevoCondition_Temperature20_Replicate$(rep)_cycle$(cycle).csv"),
            DataFrame)
        t = collect(skipmissing(df[:, 1])) ./ 24
        x = collect(skipmissing(df[:, 2]))
        if data == "host"
            push!(t_H, t); push!(H, x)
        else
            push!(t_V, t); push!(V, x)
        end
    end
end

# Indices par cycle: cycle c → reps (c-1)*3+1 : c*3
cycle_rep_idx(c) = (c-1)*3+1 : c*3

# t_all_rep par réplicat (union des temps H et V)
t_all_rep = [sort(unique(vcat(t_H[i], t_V[i]))) for i in eachindex(t_H)]

# t_end de chaque cycle = dernier temps observé parmi les 3 réplicats
t_end_cycle = [maximum(t_all_rep[i][end] for i in cycle_rep_idx(c)) for c in cycles]

# t_all_sorted par cycle (union des temps sur les 3 réplicats du cycle)
t_all_sorted_cycle = [sort(unique(vcat([t_all_rep[i] for i in cycle_rep_idx(c)]...))) for c in cycles]

# Y0 du cycle 1
function mean_first(vec_of_vecs, idxs)
    mean(first(vec_of_vecs[i]) for i in idxs)
end

Y0_global = log.([mean_first(H, cycle_rep_idx(1)), 1e-6, 1e-6, mean_first(V, cycle_rep_idx(1)), 1e-6, 1e-6, 1e-6])

# Moyennes expérimentales H0 et V0 au début de chaque cycle (pour la dilution)
H0_mean = [mean_first(H, cycle_rep_idx(c)) for c in cycles]
V0_mean = [mean_first(V, cycle_rep_idx(c)) for c in cycles]

# --- Model import ---
include(joinpath(@__DIR__, "models.jl"))
model1 = MODELS["SIVi_mu"]
model2 = MODELS["SIVi_mu_2"]


## ===== Constants =====
const K = 9.784708604680645e7


## ===== KNOWN PARAMS =====
const KNOWN_PARAMS       = [0.538001865245161, 1.0814785163571014e-8, 69.99999999999996, 0.00010000000000000026, 39.99999999999998]
const KNOWN_PARAMS_OTHERS = [1e-100]


## ===== Utils =====
function reconstruct(
    logtheta :: AbstractVector{<:Real},
    model    :: ModelSpec,
    fixed    :: Dict{Symbol,Float64} = Dict{Symbol,Float64}()
)
    fitted = Dict{Symbol,Float64}()
    idx = 1
    for p in model.fit_params
        if haskey(fixed, p)
            fitted[p] = fixed[p]
        else
            fitted[p] = exp(logtheta[idx])
            idx += 1
        end
    end
    return Float64[
        p == :k          ? K        :
        haskey(fixed, p) ? fixed[p] :
        get(fitted, p, 0.0)
        for p in model.full_params
    ]
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

function compute_error(sol, j, t_cut, t_H, H, t_V, V)
    err = 0.0
    idx_H = findall(t -> t <= t_cut, t_H[j])
    if !isempty(idx_H)
        H_pred = [logsumexp(sol(t)[1:3]) for t in t_H[j][idx_H]]
        err += sum((H_pred .- log.(H[j][idx_H])).^2)
    end
    idx_V = findall(t -> t <= t_cut, t_V[j])
    if !isempty(idx_V)
        V_pred = [logsumexp(sol(t)[4:6]) for t in t_V[j][idx_V]]
        err += sum((V_pred .- log.(V[j][idx_V])).^2)
    end
    return err
end

function compute_residue(sol, j, t1_loc, t2_loc, t_H, H, t_V, V)
    err = 0.0
    idx_H = findall(t -> t > t1_loc && t <= t2_loc, t_H[j])
    if !isempty(idx_H)
        H_pred = [logsumexp(sol(t)[1:3]) for t in t_H[j][idx_H]]
        err += sum((H_pred .- log.(H[j][idx_H])).^2)
    end
    idx_V = findall(t -> t > t1_loc && t <= t2_loc, t_V[j])
    if !isempty(idx_V)
        V_pred = [logsumexp(sol(t)[4:6]) for t in t_V[j][idx_V]]
        err += sum((V_pred .- log.(V[j][idx_V])).^2)
    end
    return err
end

function compute_derivative(t, y)
    n = length(y)
    dy = zeros(n-1)
    for i in 1:n-1
        dt = t[i+1] - t[i]
        dy[i] = (y[i+1] - y[i]) / dt
    end
    return dy
end

# --- Colors ---
host_palette  = cgrad([:darkgreen, :chartreuse])
virus_palette = cgrad([:darkred, :orangered])
col_S   = host_palette[0.0];  col_I   = host_palette[0.5];  col_R    = host_palette[1.0]
col_Vi  = virus_palette[0.0]; col_Vdp = virus_palette[0.5]; col_Vdip = virus_palette[1.0]
col_Ev  = :blue

# Dilution logic: compute Y0 for cycle c+1 from sol_end of cycle c and experimental means
function dilution_Y0(sol_end_u, H0_next, V0_next)
    S_end  = exp(sol_end_u[1]); I_end  = exp(sol_end_u[2]); R_end   = exp(sol_end_u[3])
    Vi_end = exp(sol_end_u[4]); Vdp_end= exp(sol_end_u[5]); Vdip_end= exp(sol_end_u[6])

    S_active = S_end >= 2e-6; I_active = I_end >= 2e-6; R_active = R_end >= 2e-6
    if !S_active && !I_active && !R_active
        prop_S, prop_I, prop_R = 1.0, 0.0, 0.0
    else
        H_active = (S_active ? S_end : 0.0) + (I_active ? I_end : 0.0) + (R_active ? R_end : 0.0)
        prop_S = S_active ? S_end / H_active : 0.0
        prop_I = I_active ? I_end / H_active : 0.0
        prop_R = R_active ? R_end / H_active : 0.0
    end

    Vi_active = Vi_end >= 2e-6; Vdp_active = Vdp_end >= 2e-6; Vdip_active = Vdip_end >= 2e-6
    if !Vi_active && !Vdp_active && !Vdip_active
        prop_Vi, prop_Vdp, prop_Vdip = 1.0, 0.0, 0.0
    else
        V_active = (Vi_active ? Vi_end : 0.0) + (Vdp_active ? Vdp_end : 0.0) + (Vdip_active ? Vdip_end : 0.0)
        prop_Vi   = Vi_active   ? Vi_end   / V_active : 0.0
        prop_Vdp  = Vdp_active  ? Vdp_end  / V_active : 0.0
        prop_Vdip = Vdip_active ? Vdip_end / V_active : 0.0
    end

    S0_next    = prop_S   > 0.0 ? prop_S   * H0_next : 1e-6
    I0_next    = prop_I   > 0.0 ? prop_I   * H0_next : 1e-6
    R0_next    = prop_R   > 0.0 ? prop_R   * H0_next : 1e-6
    Vi0_next   = prop_Vi  > 0.0 ? prop_Vi  * V0_next : 1e-6
    Vdp0_next  = prop_Vdp > 0.0 ? prop_Vdp * V0_next : 1e-6
    Vdip0_next = prop_Vdip> 0.0 ? prop_Vdip* V0_next : 1e-6

    return log.([S0_next, I0_next, R0_next, Vi0_next, Vdp0_next, Vdip0_next, 1e-6])
end


## ===== FitResult struct =====
struct FitResult
    cycle     :: Int
    fit_id    :: Int
    t_start   :: Float64
    t_end     :: Float64
    tvec      :: Vector{Float64}
    S    :: Vector{Float64}
    I    :: Vector{Float64}
    R    :: Vector{Float64}
    Vi   :: Vector{Float64}
    Vdp  :: Vector{Float64}
    Vdip :: Vector{Float64}
    Ev   :: Vector{Float64}
    H    :: Vector{Float64}
    V    :: Vector{Float64}
    t_residue :: Vector{Float64}
    residue   :: Vector{Float64}
    dresidue  :: Vector{Float64}
    idx_cut   :: Int
    theta_opt :: Vector{Float64}
end


## ===== Threshold / phi lists per cycle =====

# Each list has one entry per iterative model2 segment within the cycle.
# phi0_list[cyc][k] : point de départ (valeur réelle, pas log) du k-ième fit model2
# lb_list[cyc][k]   : borne inférieure (valeur réelle)
# ub_list[cyc][k]   : borne supérieure (valeur réelle)
# threshold_list[cyc][k] : seuil de résidu pour couper le segment k

threshold_list = [
    [20, 5, 4000, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20],  # cycle 1
    [1, 2, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20],  # cycle 2
    [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20],  # cycle 3
    [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20],  # cycle 4
    [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20],  # cycle 5
]

phi0_list = [
    [0, 8e-9, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20],  # cycle 1
    [0, 1e-8, 1e-8, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20, 1e-20],  # cycle 2
    [0, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9],  # cycle 3
    [0, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9],  # cycle 4
    [0, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9, 9.5e-9],  # cycle 5
]

lb_list = [
    [0, 7e-9, 1e-25, 1e-25, 1e-25, 1e-25, 1e-25, 1e-25, 1e-25, 1e-25, 1e-25, 1e-25],  # cycle 1
    [0, 1e-25, 1e-9, 1e-25, 1e-25, 1e-25, 1e-25, 1e-25, 1e-25, 1e-25, 1e-25, 1e-25],  # cycle 2
    [0, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9],  # cycle 3
    [0, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9],  # cycle 4
    [0, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9, 9e-9],  # cycle 5
]

ub_list = [
    [0, 9e-9, 1e-18, 1e-18, 1e-18, 1e-18, 1e-18, 1e-18, 1e-18, 1e-18, 1e-18, 1e-18],  # cycle 1
    [0, 1e-5, 1e-6, 1e-18, 1e-18, 1e-18, 1e-18, 1e-18, 1e-18, 1e-18, 1e-18, 1e-18],  # cycle 2
    [0, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8],  # cycle 3
    [0, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8],  # cycle 4
    [0, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8],  # cycle 5
]


## ===== Fixed params (from fit 1) =====
theta_opt1 = KNOWN_PARAMS
theta_fixed = Dict{Symbol,Float64}(
    :μ => theta_opt1[1],
    :β => theta_opt1[3],
    :δ => theta_opt1[4],
    :η => theta_opt1[5]
)


## ===== CYCLE LOOP =====

all_cycle_results = Vector{Vector{FitResult}}()
phi_segments = NamedTuple{(:t_start, :t_end, :phi), Tuple{Float64,Float64,Float64}}[]
global sol_end_u_prev = nothing
global last_phi = phi0_list[1][1]

for cyc in cycles
    println("\n" * "="^60)
    println("CYCLE $cyc")
    println("="^60)

    rep_idxs     = cycle_rep_idx(cyc)
    t_sorted     = t_all_sorted_cycle[cyc]
    t_cyc_end    = t_end_cycle[cyc]
    t_cyc_start  = t_sorted[1]
    fit_id_iter  = 1

    # --- Y0 for this cycle ---
    if cyc == 1
        Y0_cyc = Y0_global
    else
        Y0_cyc = dilution_Y0(sol_end_u_prev, H0_mean[cyc], V0_mean[cyc])
    end

    cycle_results = FitResult[]

    # ---- CYCLE 1: fit model1 first, then iterative model2 ----
    if cyc == 1
        # --- Fit 1 (model1) ---
        theta0_fit1 = log.([0.5,   1.3e-7, 50.0, 1e-3, 10.0])
        lb_fit1     = log.([0.4,   1e-7,   45.0, 1e-5,  1.0])
        ub_fit1     = log.([0.7,   1e-6,   70.0, 1e-2, 40.0])

        function obj_fit1(theta)
            p   = reconstruct(theta, model1)
            err = 0.0
            for i in rep_idxs
                prob = ODEProblem(model1.dynamics!, Y0_cyc, (t_cyc_start, t_cyc_end), reconstruct(theta, model1))
                sol  = solve(prob, Rodas5())
                err += compute_error(sol, i, t_cyc_end, t_H, H, t_V, V)
            end
            return err
        end

        if RUN_FIT_1
            println("Fit 1 - optimisation en cours...")
            res_fit1 = optimize(obj_fit1, lb_fit1, ub_fit1, theta0_fit1, Fminbox(BFGS()))
            theta_opt1_loc = exp.(Optim.minimizer(res_fit1))
        else
            println("Fit 1 - parametres charges manuellement.")
            theta_opt1_loc = KNOWN_PARAMS
        end
        println("mu=$(theta_opt1_loc[1]), phi=$(theta_opt1_loc[2]), beta=$(theta_opt1_loc[3]), delta=$(theta_opt1_loc[4]), eta=$(theta_opt1_loc[5])")

        # Simulate fit1
        p_fit1    = reconstruct(log.(theta_opt1_loc), model1)
        t1_dense  = range(t_cyc_start, t_cyc_end, length=500)
        prob_fit1 = ODEProblem(model1.dynamics!, Y0_cyc, (t_cyc_start, t_cyc_end), p_fit1)
        sol_fit1  = solve(prob_fit1, Rodas5(); saveat=t1_dense)
        S1,I1,R1,Vi1,Vdp1,Vdip1,Ev1 = extract_all(sol_fit1, t1_dense)
        H1 = S1 .+ I1 .+ R1
        V1 = Vi1 .+ Vdp1 .+ Vdip1

        # Residue for fit1
        residue1 = Float64[]
        for i in 2:length(t_sorted)
            err = sum(compute_residue(sol_fit1, j, 0.0, t_sorted[i], t_H, H, t_V, V) for j in rep_idxs)
            push!(residue1, err)
        end
        i_threshold1 = findlast(<(threshold_list[1][1]), residue1)
        t_threshold1 = t_sorted[i_threshold1 + 1]
        println("Fit 1 threshold = $(threshold_list[1][1]) => t_threshold = $t_threshold1")

        mask1     = t1_dense .<= t_threshold1
        idx_end1  = findlast(mask1)
        t1_vec    = collect(t1_dense[mask1])

        # Store fit1 as FitResult
        t_residue1_vec = t_sorted[2:end]
        dresidue1      = compute_derivative(t_residue1_vec, residue1)
        push!(cycle_results, FitResult(
            1, 1,
            t_cyc_start, t_threshold1,
            t1_vec,
            S1[mask1], I1[mask1], R1[mask1],
            Vi1[mask1], Vdp1[mask1], Vdip1[mask1], Ev1[mask1],
            H1[mask1], V1[mask1],
            t_residue1_vec, residue1, dresidue1,
            i_threshold1,
            theta_opt1_loc
        ))

        push!(phi_segments, (t_start=t_cyc_start, t_end=t_threshold1, phi=theta_opt1_loc[2]))
        global Y0_iter    = log.([S1[idx_end1], I1[idx_end1], R1[idx_end1], Vi1[idx_end1], Vdp1[idx_end1], Vdip1[idx_end1], Ev1[idx_end1]])
        global t_start_iter = t_threshold1
        global last_phi = phi0_list[1][1]
        global fit_id_iter += 1

    else
        # For cycles 2-5: start directly with model2
        global Y0_iter      = Y0_cyc
        global t_start_iter = t_cyc_start
        global fit_id_iter += 1
    end

    # ---- Iterative model2 fits ----

    while t_start_iter < t_cyc_end && fit_id_iter <= length(threshold_list[cyc])
        println("\n--- Cycle $cyc | Fit $fit_id_iter  t_start = $t_start_iter ---")

        threshold = threshold_list[cyc][fit_id_iter]
        phi0_fi   = (fit_id_iter == 1) ? last_phi : phi0_list[cyc][fit_id_iter]
        lb_fi_v   = lb_list[cyc][fit_id_iter]
        ub_fi_v   = ub_list[cyc][fit_id_iter]

        # sécurité : clamp theta0 dans les bornes
        phi0_fi = clamp(phi0_fi, lb_fi_v, ub_fi_v)

        if RUN_FIT_OTHERS
            # Si on est dans un cycle > 1 et qu'on gèle phi : pas d'optimisation
            if !FIT_PHI_NEXT_CYCLES && cyc > 1
                theta_opt_fi = [last_phi]   # on réutilise le dernier phi du cycle 1
            else
                theta0_fi_log = [log(phi0_fi)]
                lb_fi         = [log(lb_fi_v)]
                ub_fi         = [log(ub_fi_v)]

                Y0_loc      = copy(Y0_iter)
                t_start_loc = t_start_iter

                obj_fi = function(theta)
                    p    = reconstruct(theta, model2, theta_fixed)
                    prob = ODEProblem(model2.dynamics!, Y0_loc, (t_start_loc, t_cyc_end), p)
                    sol  = solve(prob, Rodas5())
                    err  = 0.0
                    for j in rep_idxs
                        err += compute_error(sol, j, t_cyc_end, t_H, H, t_V, V)
                    end
                    println(err)
                    return err
                end

                res_fi       = optimize(obj_fi, lb_fi, ub_fi, theta0_fi_log, Fminbox(BFGS()))
                println(Optim.converged(res_fi))
                println(Optim.minimum(res_fi))
                println(Optim.iterations(res_fi))
                theta_opt_fi = exp.(Optim.minimizer(res_fi))
            end
        else
            theta_opt_fi = KNOWN_PARAMS_OTHERS
        end

        println("  phi_i = $(theta_opt_fi[1])")
        global last_phi = theta_opt_fi[1]

        # Simulate
        t_fi_dense = range(t_start_iter, t_cyc_end, length=500)
        p_fi       = reconstruct(log.(theta_opt_fi), model2, theta_fixed)
        prob_fi    = ODEProblem(model2.dynamics!, Y0_iter, (t_start_iter, t_cyc_end), p_fi)
        sol_fi     = solve(prob_fi, Rodas5(); saveat=t_fi_dense)

        S_fi,I_fi,R_fi,Vi_fi,Vdp_fi,Vdip_fi,Ev_fi = extract_all(sol_fi, t_fi_dense)
        H_fi = S_fi .+ I_fi .+ R_fi
        V_fi = Vi_fi .+ Vdp_fi .+ Vdip_fi

        # Residue
        idx_start_loc = findfirst(t -> t >= t_start_iter, t_sorted)
        t_fi_data = isnothing(idx_start_loc) || idx_start_loc == 1 ? t_sorted : t_sorted[idx_start_loc-1:end]

        residue_local = Float64[]
        for i in 2:length(t_fi_data)
            err = sum(compute_residue(sol_fi, j, t_fi_data[i-1], t_fi_data[i], t_H, H, t_V, V) for j in rep_idxs)
            push!(residue_local, err)
        end
        t_residue_local = t_fi_data[2:end]
        dresidue        = compute_derivative(t_residue_local, residue_local)

        # Threshold detection
        idx_local = findlast(x -> x < threshold, residue_local)
        t_end_cur = isnothing(idx_local) ? t_sorted[end] : t_fi_data[idx_local + 1]
        println("  t_end = $t_end_cur")

        # Extract segment
        idx_end_dense = something(findfirst(t -> t >= t_end_cur, t_fi_dense), length(t_fi_dense))
        t_plot = collect(t_fi_dense[1:idx_end_dense])
        n      = length(t_plot)

        idx_cut = isnothing(idx_local) ? length(residue_local) : idx_local

        push!(cycle_results, FitResult(
            cyc, fit_id_iter,
            t_start_iter, t_end_cur,
            t_plot,
            S_fi[1:n], I_fi[1:n], R_fi[1:n],
            Vi_fi[1:n], Vdp_fi[1:n], Vdip_fi[1:n], Ev_fi[1:n],
            H_fi[1:n], V_fi[1:n],
            t_residue_local, residue_local, dresidue,
            idx_cut,
            theta_opt_fi
        ))

        push!(phi_segments, (t_start=t_start_iter, t_end=t_end_cur, phi=theta_opt_fi[1]))

        # Update
        Y0_iter = log.([S_fi[idx_end_dense], I_fi[idx_end_dense], R_fi[idx_end_dense], Vi_fi[idx_end_dense], Vdp_fi[idx_end_dense], Vdip_fi[idx_end_dense], Ev_fi[idx_end_dense]])
        t_start_iter = t_end_cur
        fit_id_iter += 1

        t_start_iter >= t_cyc_end && break
    end

    # Save sol end for dilution
    # Reconstruct sol at t_cyc_end using last fit
    last_fr = cycle_results[end]
    p_last  = cyc == 1 && length(cycle_results) == 1 ?
              reconstruct(log.(last_fr.theta_opt), model1) :
              reconstruct(log.(last_fr.theta_opt), model2, theta_fixed)
    model_last = (cyc == 1 && length(cycle_results) == 1) ? model1 : model2
    prob_last  = ODEProblem(model_last.dynamics!, log.(vcat(
                    [last_fr.S[1], last_fr.I[1], last_fr.R[1],
                     last_fr.Vi[1], last_fr.Vdp[1], last_fr.Vdip[1], last_fr.Ev[1]])),
                    (last_fr.t_start, t_cyc_end), p_last)
    sol_last       = solve(prob_last, Rodas5())
    global sol_end_u_prev = sol_last.u[end]

    push!(all_cycle_results, cycle_results)
    println("\nCycle $cyc done: $(length(cycle_results)) fit segments")
end


## ===== PLOTTING =====

run_id   = Dates.format(t_global_start, "yyyymmdd-HHMMSS")
if FIT_PHI_NEXT_CYCLES
    out_dir  = joinpath(@__DIR__, "140426_output")
else
    out_dir  = joinpath(@__DIR__, "140426_output_keep_phi")
end
mkpath(out_dir)

# Global y limits across all data and simulations
all_vals_global = vcat(
    vcat(H...), vcat(V...),
    [vcat(fr.H, fr.V, fr.S, fr.I, fr.R, fr.Vi, fr.Vdp, fr.Vdip, fr.Ev)
     for cres in all_cycle_results for fr in cres]...
)
y_min_g = max(1e-5, minimum(filter(x -> isfinite(x) && x > 0, all_vals_global)))
y_max_g = maximum(filter(x -> isfinite(x) && x > 0, all_vals_global))

# ── Per-cycle plot (one row per cycle, 3 cols: H+V, components, residue) ──
n_rows_cyc = 5
fig_cyc = plot(layout=(n_rows_cyc, 3), grid=true,
               size=(2800, 800 * n_rows_cyc),
               xlabel="Time (days)", legend=:bottomleft, margins=15mm)

for cyc in cycles
    rep_idxs  = cycle_rep_idx(cyc)
    cres      = all_cycle_results[cyc]
    t_cyc_end = t_end_cycle[cyc]
    row       = cyc

    col_hv  = (row-1)*3 + 1
    col_cmp = (row-1)*3 + 2
    col_res = (row-1)*3 + 3

    # Scatter: only data points of this cycle
    for i in rep_idxs
        scatter!(fig_cyc[col_hv], t_H[i], H[i], color=:green, marker=:circle, alpha=0.5,
                 label=(i==rep_idxs[1] ? "H data" : false))
        scatter!(fig_cyc[col_hv], t_V[i], V[i], color=:red,   marker=:square, alpha=0.5,
                 label=(i==rep_idxs[1] ? "V data" : false))
    end

    for fr in cres
        lab = "fit $(fr.fit_id)"
        plot!(fig_cyc[col_hv], fr.tvec, fr.H, lw=3, color=:green, label="H $lab", yscale=:log10)
        plot!(fig_cyc[col_hv], fr.tvec, fr.V, lw=3, color=:red,   label="V $lab")

        plot!(fig_cyc[col_cmp], fr.tvec, fr.S,    lw=2, label="S $lab",    color=col_S,    yscale=:log10)
        plot!(fig_cyc[col_cmp], fr.tvec, fr.I,    lw=2, label="I $lab",    color=col_I)
        plot!(fig_cyc[col_cmp], fr.tvec, fr.R,    lw=2, label="R $lab",    color=col_R)
        plot!(fig_cyc[col_cmp], fr.tvec, fr.Vi,   lw=2, label="Vi $lab",   color=col_Vi)
        plot!(fig_cyc[col_cmp], fr.tvec, fr.Vdp,  lw=2, label="Vdp $lab",  color=col_Vdp)
        plot!(fig_cyc[col_cmp], fr.tvec, fr.Vdip, lw=2, label="Vdip $lab", color=col_Vdip)
        plot!(fig_cyc[col_cmp], fr.tvec, fr.Ev,   lw=2, label="Ev $lab",   color=col_Ev)

        plot!(fig_cyc[col_res], fr.t_residue[1:fr.idx_cut], fr.residue[1:fr.idx_cut], lw=2, label="residue $lab", color=:black)
        #if length(fr.dresidue) > 0
        #    plot!(fig_cyc[col_res], fr.t_residue[1:fr.idx_cut-1], fr.dresidue[1:fr.idx_cut-1], lw=2, ls=:dash,
        #          label="d res/dt $lab", color=:grey)
        #end
        vline!(fig_cyc[col_res], [fr.t_end], lw=1, ls=:dot, color=:blue, label=false)
    end

    # Vertical separators between fit segments
    for fr in cres[1:end-1]
        vline!(fig_cyc[col_hv],  [fr.t_end], lw=1, ls=:dot, color=:black, label=false)
        vline!(fig_cyc[col_cmp], [fr.t_end], lw=1, ls=:dot, color=:black, label=false)
    end

    ylabel!(fig_cyc[col_hv],  "Concentration (parts/mL)")
    ylabel!(fig_cyc[col_cmp], "Concentration (parts/mL)")
    ylabel!(fig_cyc[col_res], "Residue")
    title!(fig_cyc[col_hv],  "Cycle $cyc - H and V")
    title!(fig_cyc[col_cmp], "Cycle $cyc - S,I,R,Vi,Vdp,Vdip,Ev")
    title!(fig_cyc[col_res], "Cycle $cyc - Residue")
    ylims!(fig_cyc[col_hv],  (y_min_g, y_max_g))
    ylims!(fig_cyc[col_cmp], (y_min_g, y_max_g))
    xlims!(fig_cyc[col_hv],  (t_end_cycle[1] * (cyc-1) / 5, t_cyc_end))  # rough x range
    xlims!(fig_cyc[col_cmp], (t_end_cycle[1] * (cyc-1) / 5, t_cyc_end))
    xlims!(fig_cyc[col_res], (t_end_cycle[1] * (cyc-1) / 5, t_cyc_end))
end

cyc_plot_path = joinpath(out_dir, "$(run_id)_per_cycle_plot.png")
savefig(fig_cyc, cyc_plot_path)
println("Per-cycle plot saved: $cyc_plot_path")

# ── All cycles combined plot (1 row, 3 cols) ──
fig_all = plot(layout=(1, 3), grid=true, size=(2200, 800),
               xlabel="Time (days)", legend=:bottomleft, margins=15mm)

for cyc in cycles
    rep_idxs = cycle_rep_idx(cyc)
    cres     = all_cycle_results[cyc]

    for i in rep_idxs
        scatter!(fig_all[1], t_H[i], H[i], color=:green, marker=:circle, alpha=0.4,
                 label=(cyc==1 && i==rep_idxs[1] ? "H data" : false))
        scatter!(fig_all[1], t_V[i], V[i], color=:red,   marker=:square, alpha=0.4,
                 label=(cyc==1 && i==rep_idxs[1] ? "V data" : false))
    end

    for fr in cres
        plot!(fig_all[1], fr.tvec, fr.H, lw=2, color=:green, label=false, yscale=:log10)
        plot!(fig_all[1], fr.tvec, fr.V, lw=2, color=:red,   label=false)
        plot!(fig_all[2], fr.tvec, fr.S,    lw=2, color=col_S,    label=false, yscale=:log10)
        plot!(fig_all[2], fr.tvec, fr.I,    lw=2, color=col_I,    label=false)
        plot!(fig_all[2], fr.tvec, fr.R,    lw=2, color=col_R,    label=false)
        plot!(fig_all[2], fr.tvec, fr.Vi,   lw=2, color=col_Vi,   label=false)
        plot!(fig_all[2], fr.tvec, fr.Vdp,  lw=2, color=col_Vdp,  label=false)
        plot!(fig_all[2], fr.tvec, fr.Vdip, lw=2, color=col_Vdip, label=false)
        plot!(fig_all[2], fr.tvec, fr.Ev,   lw=2, color=col_Ev,   label=false)
        plot!(fig_all[3], fr.t_residue[1:fr.idx_cut], fr.residue[1:fr.idx_cut], lw=2, color=:black, label=false)
        #plot!(fig_all[3], fr.t_residue[1:fr.idx_cut-1], fr.dresidue[1:fr.idx_cut-1], lw=2, color=:black, label=false, ls=:dash)
    end

    # Mark cycle boundaries
    if cyc < 5
        vline!(fig_all[1], [t_end_cycle[cyc]], lw=1.5, ls=:dash, color=:navy,  label=false)
        vline!(fig_all[2], [t_end_cycle[cyc]], lw=1.5, ls=:dash, color=:navy,  label=false)
        vline!(fig_all[3], [t_end_cycle[cyc]], lw=1.5, ls=:dash, color=:navy,  label=false)
    end
end

for seg in phi_segments
    if !(seg.t_end in t_end_cycle)
        vline!(fig_all[1], [seg.t_end], lw=1, ls=:dot, color=:black, label=false)
        vline!(fig_all[2], [seg.t_end], lw=1, ls=:dot, color=:black, label=false)
        vline!(fig_all[3], [seg.t_end], lw=1, ls=:dot, color=:black, label=false)
    end
end

# Legend entries for components
plot!(fig_all[2], [], [], lw=2, color=col_S,    label="S")
plot!(fig_all[2], [], [], lw=2, color=col_I,    label="I")
plot!(fig_all[2], [], [], lw=2, color=col_R,    label="R")
plot!(fig_all[2], [], [], lw=2, color=col_Vi,   label="Vi")
plot!(fig_all[2], [], [], lw=2, color=col_Vdp,  label="Vdp")
plot!(fig_all[2], [], [], lw=2, color=col_Vdip, label="Vdip")
plot!(fig_all[2], [], [], lw=2, color=col_Ev,   label="Ev")

ylabel!(fig_all[1], "Concentration (parts/mL)")
ylabel!(fig_all[2], "Concentration (parts/mL)")
ylabel!(fig_all[3], "Residue")
title!(fig_all[1], "All cycles - H and V")
title!(fig_all[2], "All cycles - S,I,R,Vi,...")
title!(fig_all[3], "All cycles - Residue")
ylims!(fig_all[1], (y_min_g, y_max_g))
ylims!(fig_all[2], (y_min_g, y_max_g))

all_plot_path = joinpath(out_dir, "$(run_id)_all_cycles_plot.png")
savefig(fig_all, all_plot_path)
println("All-cycles plot saved: $all_plot_path")

# ── phi vs time (piecewise constant) ──
fig_phi = plot(grid=true, size=(1200, 400),
               xlabel="Time (days)", ylabel="phi",
               title="phi vs time (piecewise constant)",
               legend=false, margins=15mm)

for seg in phi_segments
    plot!(fig_phi, [seg.t_start, seg.t_end], [seg.phi, seg.phi], lw=3, color=:purple, yscale=:log10)
    # vertical connector to next segment value (step)
end
# Add step connectors between consecutive segments
for k in 1:length(phi_segments)-1
    t_conn = phi_segments[k].t_end
    phi_a  = phi_segments[k].phi
    phi_b  = phi_segments[k+1].phi
    plot!(fig_phi, [t_conn, t_conn], [phi_a, phi_b], lw=1, ls=:dot, color=:grey)
end
# Cycle boundary markers
for cyc in cycles
    vline!(fig_phi, [t_end_cycle[cyc]], lw=1.5, ls=:dash, color=:navy, label=false)
end

phi_plot_path = joinpath(out_dir, "$(run_id)_phi_vs_time.png")
savefig(fig_phi, phi_plot_path)
println("Phi plot saved: $phi_plot_path")


## ===== Save log =====
log_path = joinpath(out_dir, "$(run_id)_log.txt")

open(log_path, "w") do io
    sep = "=" ^ 60
    println(io, sep)
    println(io, "  RUN LOG - iterative fits over 5 cycles")
    println(io, "  Date      : $(Dates.format(t_global_start, "yyyy-mm-dd HH:MM:SS"))")
    println(io, "  RUN_FIT_1 : $(RUN_FIT_1 ? "optimisation" : "parametres manuels")")
    println(io, "  RUN_FIT_OTHERS : $RUN_FIT_OTHERS")
    println(io, "  FIT_PHI_NEXT_CYCLES : $FIT_PHI_NEXT_CYCLES")
    println(io, sep)

    println(io, "\n--- FIXED PARAMS (from fit 1) ---")
    for (p, v) in theta_fixed
        println(io, "  $(rpad(string(p), 8)) = $(round(v, sigdigits=6))")
    end

    for cyc in cycles
        println(io, "\n$sep")
        println(io, "  CYCLE $cyc  (t_end = $(t_end_cycle[cyc]) d)")
        println(io, sep)
        cres = all_cycle_results[cyc]
        for fr in cres
            println(io, "\n  Fit $(fr.fit_id)")
            println(io, "    Segment  : $(round(fr.t_start, digits=4)) to $(round(fr.t_end, digits=4)) d")
            if length(fr.theta_opt) == 1
                println(io, "    phi_i    = $(round(fr.theta_opt[1], sigdigits=6))")
            else
                for (k, p) in enumerate(model1.fit_params)
                    println(io, "    $(rpad(string(p), 8)) = $(round(fr.theta_opt[k], sigdigits=6))")
                end
            end
            println(threshold_list)
            println(io, "    threshold = $(threshold_list[cyc][fr.fit_id])")
        end
    end

    println(io, "\n$sep")
    println(io, "  PHI SEGMENTS")
    println(io, sep)
    for seg in phi_segments
        println(io, "  t=[$(round(seg.t_start,digits=3)), $(round(seg.t_end,digits=3))]  phi=$(round(seg.phi, sigdigits=6))")
    end

    t_total = (now() - t_global_start).value / 1000
    println(io, "\n$sep")
    println(io, "  Total runtime : $(round(t_total, digits=1))s")
    println(io, "  END OF LOG")
    println(io, sep)
end

println("\n" * "="^60)
println("ALL DONE in $(round((now() - t_global_start).value / 1000, digits=1))s")
println("Log:           $log_path")
println("Per-cycle:     $cyc_plot_path")
println("All cycles:    $all_plot_path")
println("Phi vs time:   $phi_plot_path")
println("="^60 * "\n")