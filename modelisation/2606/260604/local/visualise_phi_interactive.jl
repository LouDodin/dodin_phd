using CSV
using DataFrames
using DifferentialEquations
using SciMLBase
using OrdinaryDiffEqRosenbrock
using DataInterpolations
using Statistics
using GLMakie
using Printf

# ─────────────────────────────────────────────────────────────────────────────
# plot_interactive.jl
#
# Interactive visualisation with one slider per cycle (n₁ … n₅ interior nodes).
# Turning a slider instantly reloads the corresponding polynomial.txt and
# redraws φ(t), S(t) and V(t) on the three panels.
#
# Expected file layout (same as plot_fit.jl):
#   <output_root>/nint_<n1>_<n2>_<n3>_<n4>_<n5>/polynomial.txt
#
# Usage:
#   julia plot_interactive.jl
#   julia plot_interactive.jl /path/to/output_root /path/to/data_dir
# ─────────────────────────────────────────────────────────────────────────────

# ── Paths ─────────────────────────────────────────────────────────────────────
const OUTPUT_ROOT = length(ARGS) >= 1 ? ARGS[1] :
                    joinpath(@__DIR__, "../genotoul/output")
const DATA_DIR    = length(ARGS) >= 2 ? ARGS[2] :
                    joinpath(@__DIR__, "../../../input/xp_input_20")

# ── Biology ───────────────────────────────────────────────────────────────────
const r = 0.574619342477644
const K = 6.675449070379925e7
const β = 144.0
const δ = 0.02

const REPLICATES = ["A", "B", "C"]
const N_CYCLES   = 5
const NODE_RANGE = 0:5          # allowed values per slider

# ── Data structures ───────────────────────────────────────────────────────────
struct Cycle
    rep   :: String
    index :: Int
    tH    :: Vector{Float64}
    H     :: Vector{Float64}
    tV    :: Vector{Float64}
    V     :: Vector{Float64}
    u0    :: Vector{Float64}
end

# ── Load experimental data once ───────────────────────────────────────────────
function load_replicate(rep::String)::Vector{Cycle}
    cycles     = Cycle[]
    t_offset   = nothing
    t_end_prev = nothing
    for cyc in 1:N_CYCLES
        path_H = joinpath(DATA_DIR,
            "hostData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cyc).csv")
        path_V = joinpath(DATA_DIR,
            "virusData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cyc).csv")
        (!isfile(path_H) || !isfile(path_V)) && continue
        df_H = CSV.read(path_H, DataFrame)
        df_V = CSV.read(path_V, DataFrame)
        tH = df_H[:, 1] ./ 24.0;  H = df_H[:, 2]
        tV = df_V[:, 1] ./ 24.0;  V = df_V[:, 2]
        if t_offset === nothing; t_offset = tH[1]; end
        tH .-= t_offset;  tV .-= t_offset
        if t_end_prev !== nothing
            gap = tH[1] - t_end_prev
            tH .-= gap;  tV .-= gap
        end
        t_end_prev = tH[end]
        push!(cycles, Cycle(rep, cyc, tH, H, tV, V, [H[1], V[1]]))
    end
    return cycles
end

println("Loading experimental data...")
all_cycles     = vcat([load_replicate(rep) for rep in REPLICATES]...)
cycles_per_rep = Dict(rep => filter(c -> c.rep == rep, all_cycles) for rep in REPLICATES)
println("$(length(all_cycles)) cycles loaded.")

t_global_min = minimum(c.tH[1]   for c in all_cycles)
t_global_max = maximum(c.tH[end] for c in all_cycles)

# ── Parse polynomial.txt ──────────────────────────────────────────────────────
function parse_knots(path::String)
    lines = readlines(path)
    t_knots = Float64[];  θ = Float64[]
    in_knots = false
    for line in lines
        if occursin("OPTIMAL KNOTS", line);  in_knots = true;  continue;  end
        in_knots || continue
        isempty(strip(line)) && break
        occursin(r"^\s*i\s", line) && continue
        parts = split(strip(line))
        length(parts) >= 3 || continue
        try
            push!(t_knots, parse(Float64, parts[2]))
            push!(θ,       parse(Float64, parts[3]))
        catch; continue; end
    end
    isempty(t_knots) && error("No knots found in $path")
    return θ, t_knots
end

# ── ODE ───────────────────────────────────────────────────────────────────────
function ode_model!(dY, Y, p, t)
    S, Vi = Y[1], Y[2]
    ϕ = p(t)
    dY[1] = r * S * (1 - S/K) - ϕ * S * Vi
    dY[2] = β * ϕ * S * Vi - δ * Vi
end

function integrate_cycle(u0, t0, t1, phi_func; n_points=400)
    sv   = collect(range(t0, t1; length=n_points))
    prob = ODEProblem(ode_model!, u0, (t0, t1), phi_func)
    try
        sol = solve(prob, Rodas5(),
                    reltol=1e-6, abstol=1e-6, saveat=sv,
                    isoutofdomain=(u, p, t) -> any(x -> x < 0, u))
        return SciMLBase.successful_retcode(sol) ? sol : nothing
    catch
        return nothing
    end
end

# ── Compute model curves for a given combination ──────────────────────────────
# Returns (t_fine, phi_vals, sol_S_per_cycle, sol_V_per_cycle, sol_t_per_cycle)
# or nothing on failure.
function compute_combination(ns::NTuple{5,Int})
    combi = join(ns, "_")
    path  = joinpath(OUTPUT_ROOT, "nint_$(combi)", "polynomial.txt")
    if !isfile(path)
        @warn "File not found: $path"
        return nothing
    end
    θ, t_knots = parse_knots(path)
    spline   = CubicSpline(θ, t_knots)
    phi_func = t -> exp(spline(clamp(t, t_knots[1], t_knots[end])))

    t_fine   = collect(range(t_global_min, t_global_max; length=2000))
    phi_vals = phi_func.(t_fine)

    sol_S = Vector{Float64}[]
    sol_V = Vector{Float64}[]
    sol_t = Vector{Float64}[]

    for cyc_idx in 1:N_CYCLES
        cycs = filter(c -> c.index == cyc_idx, all_cycles)
        isempty(cycs) && continue
        u0_mean = [mean(c.u0[1] for c in cycs), mean(c.u0[2] for c in cycs)]
        t0 = minimum(c.tH[1]                   for c in cycs)
        t1 = maximum(max(c.tH[end], c.tV[end]) for c in cycs)
        sol = integrate_cycle(u0_mean, t0, t1, phi_func)
        if sol !== nothing
            push!(sol_t, sol.t)
            push!(sol_S, max.(sol[1, :], 1e-12))
            push!(sol_V, max.(sol[2, :], 1e-12))
        else
            push!(sol_t, Float64[])
            push!(sol_S, Float64[])
            push!(sol_V, Float64[])
        end
    end

    return (t_fine=t_fine, phi=phi_vals, t_knots=t_knots, theta=θ,
            sol_t=sol_t, sol_S=sol_S, sol_V=sol_V)
end

# ── Colours ───────────────────────────────────────────────────────────────────
REP_COLORS   = [RGBf(0.6, 0.8, 1.0), RGBf(31/255, 119/255, 180/255), RGBf(0.0, 0.3, 0.7)]
MODEL_COLOR  = RGBf(255/255, 127/255, 14/255)
PHI_COLOR    = RGBf(31/255, 119/255, 180/255)
KNOT_COLOR   = RGBf(0.9, 0.3, 0.1)

# ── Build figure ──────────────────────────────────────────────────────────────
fig = Figure(size=(1400, 1200), fontsize=16)

# Title
Label(fig[0, 1:2], "Interactive φ(t) / S(t) / V(t)";
      fontsize=18, font=:bold)

# Three plot panels (left column)
ax_phi = Axis(fig[1, 1];
    ylabel="φ(t)  [mL cell⁻¹ day⁻¹]",
    title="Infection rate φ(t)",
    yscale=log10,
    ygridvisible=true,
)

ax_S = Axis(fig[2, 1];
    ylabel="S(t)  [cells mL⁻¹]",
    title="Hosts S(t)",
    yscale=log10,
    ygridvisible=true,
)
ylims!(ax_S, 1e3, 1e9)

ax_V = Axis(fig[3, 1];
    xlabel="Time (days)",
    ylabel="V(t)  [virions mL⁻¹]",
    title="Virus V(t)",
    yscale=log10,
    ygridvisible=true,
)
ylims!(ax_V, 1e6, 1e10)


# Slider panel (right column)
slider_panel = fig[1:3, 2] = GridLayout()

Label(slider_panel[0, 1:2], "Interior nodes per cycle"; font=:bold, fontsize=15)

sliders = SliderGrid(
    slider_panel[1, 1:2],
    (label="Cycle 1  n₁", range=NODE_RANGE, startvalue=3),
    (label="Cycle 2  n₂", range=NODE_RANGE, startvalue=2),
    (label="Cycle 3  n₃", range=NODE_RANGE, startvalue=2),
    (label="Cycle 4  n₄", range=NODE_RANGE, startvalue=3),
    (label="Cycle 5  n₅", range=NODE_RANGE, startvalue=2),
    width=320,
    tellheight=false,
)

# Status label
status_label = Label(slider_panel[2, 1:2], "combination: 3_2_2_3_2";
                     fontsize=13, color=:gray40)

# Column widths: plots take ~75%, sliders ~25%
colsize!(fig.layout, 1, Relative(0.75))
colsize!(fig.layout, 2, Relative(0.25))

# ── Plot static experimental data ─────────────────────────────────────────────
# Collect one handle per replicate (first cycle) for the legends
rep_handles_S = []   # scatter handles on ax_S
rep_handles_V = []   # scatter handles on ax_V

for (i, rep) in enumerate(REPLICATES)
    col  = REP_COLORS[i]
    cycs = cycles_per_rep[rep]
    isempty(cycs) && continue
    for (j, cyc) in enumerate(cycs)
        p_S = scatter!(ax_S, cyc.tH, cyc.H;
            color=col, alpha=0.7, markersize=7, strokewidth=0)
        p_V = scatter!(ax_V, cyc.tV, cyc.V;
            color=col, alpha=0.7, markersize=7, strokewidth=0)
        if j == 1
            push!(rep_handles_S, p_S)
            push!(rep_handles_V, p_V)
        end
    end
end

# Cycle boundary lines on all three axes
for ax in (ax_phi, ax_S, ax_V)
    for cyc in cycles_per_rep["A"]
        vlines!(ax, [cyc.tH[1]]; color=:gray60, linewidth=1, linestyle=:dash)
    end
end
# Cycle labels on ax_phi
for cyc in cycles_per_rep["A"]
    text!(ax_phi, cyc.tH[1] + 0.2, 1e-8;
        text="C$(cyc.index)", fontsize=11, color=:gray50, align=(:left, :bottom))
end

# ── Observables for dynamic model curves ──────────────────────────────────────
t_phi_obs = Observable(Float64[t_global_min, t_global_max])
phi_obs   = Observable([1e-9, 1e-9])
tk_obs    = Observable(Float64[])
kv_obs    = Observable(Float64[])

sol_t_obs = [Observable(Float64[]) for _ in 1:N_CYCLES]
sol_S_obs = [Observable(Float64[]) for _ in 1:N_CYCLES]
sol_V_obs = [Observable(Float64[]) for _ in 1:N_CYCLES]

# Draw dynamic lines (φ panel)
lines!(ax_phi, t_phi_obs, phi_obs; color=PHI_COLOR, linewidth=2.5, label="φ(t)")
scatter!(ax_phi, tk_obs, kv_obs;
    color=KNOT_COLOR, markersize=8, strokewidth=0, label="knots")
axislegend(ax_phi; position=:rt, framevisible=false)

# Draw dynamic model lines — capture first handle for the legend
model_handle_S = lines!(ax_S, sol_t_obs[1], sol_S_obs[1];
    color=MODEL_COLOR, linewidth=2.5)
model_handle_V = lines!(ax_V, sol_t_obs[1], sol_V_obs[1];
    color=MODEL_COLOR, linewidth=2.5)
for cyc_idx in 2:N_CYCLES
    lines!(ax_S, sol_t_obs[cyc_idx], sol_S_obs[cyc_idx];
        color=MODEL_COLOR, linewidth=2.5)
    lines!(ax_V, sol_t_obs[cyc_idx], sol_V_obs[cyc_idx];
        color=MODEL_COLOR, linewidth=2.5)
end

# ── Legends for ax_S and ax_V ─────────────────────────────────────────────────
rep_labels = ["Rep $r" for r in REPLICATES]

axislegend(ax_S,
    [rep_handles_S..., model_handle_S],
    [rep_labels..., "Model"];
    position=:rt, framevisible=false)

axislegend(ax_V,
    [rep_handles_V..., model_handle_V],
    [rep_labels..., "Model"];
    position=:rt, framevisible=false)

# ── Callback: recompute and push to Observables ───────────────────────────────
function update_plots(ns)
    combi = join(ns, "_")
    status_label.text[] = "Computing combination: $(combi)…"
    result = compute_combination(Tuple(ns))
    if result === nothing
        status_label.text[] = "⚠  combination $(combi): file not found"
        return
    end
    # φ(t)
    t_phi_obs[] = result.t_fine
    phi_obs[]   = result.phi
    tk_obs[]    = result.t_knots
    kv_obs[]    = exp.(result.theta)
    # ODE per cycle
    for cyc_idx in 1:N_CYCLES
        if cyc_idx <= length(result.sol_t) && !isempty(result.sol_t[cyc_idx])
            sol_t_obs[cyc_idx][] = result.sol_t[cyc_idx]
            sol_S_obs[cyc_idx][] = result.sol_S[cyc_idx]
            sol_V_obs[cyc_idx][] = result.sol_V[cyc_idx]
        else
            sol_t_obs[cyc_idx][] = Float64[]
            sol_S_obs[cyc_idx][] = Float64[]
            sol_V_obs[cyc_idx][] = Float64[]
        end
    end
    status_label.text[] = "Showing combination: $(combi)  ✓"
end

# Wire sliders → callback (fires on every slider change)
slider_vals = [s.value for s in sliders.sliders]
on(events(fig).keyboardbutton) do _; end   # keep event loop alive

lift(slider_vals...) do vals...
    update_plots(collect(vals))
end

# Initial draw
update_plots([3, 2, 2, 3, 2])

# ── Launch ────────────────────────────────────────────────────────────────────
println("\nLaunching interactive window…  (close window to quit)")
display(fig)
wait(fig.scene)   # block until window is closed