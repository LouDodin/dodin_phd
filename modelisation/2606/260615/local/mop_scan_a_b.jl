## ===== Packages =====
using DifferentialEquations
using OrdinaryDiffEqRosenbrock
using Plots
using Measures
using LaTeXStrings
using Statistics
using Printf
using SciMLBase

## ===== Load model =====
const MODEL_FILE = joinpath(@__DIR__, "models/model_SRVi_scan_a_b.jl")
println("Loading model from: $MODEL_FILE")
include(MODEL_FILE)
using .ModelDef

## ===== Settings =====
const isoutofdomain = (u, p, t) -> any(x -> x < 0, u)
const prop_S_0 = 1.0
const t_span   = (0.0, 30.0)
const t_save   = 0.0:0.1:60.0
const a  = 1e-6
const b  = 1e-6

e = floor(Int, log10(a))
m = a / 10.0^e
m_str = m ≈ round(m) ? "$(Int(round(m)))" : @sprintf("%.1f", m)
a_str = "$(m_str)e$e"

e = floor(Int, log10(b))
m = b / 10.0^e
m_str = m ≈ round(m) ? "$(Int(round(m)))" : @sprintf("%.1f", m)
b_str = "$(m_str)e$e"

mop_values = 10 .^ range(-1, 10,  length=40)
H0_values  = [1e6]#[1e4, 5e4, 1e5, 5e5, 1e6, 5e6, 1e7, 5e7, 1e8]# [1e6] # [1e4, 5e4, 1e5, 5e5, 1e6, 5e6, 1e7, 5e7, 1e8]

output_path = joinpath(@__DIR__, "output/alpha vs a_b/a_b/1e6")
mkpath(output_path)

## ===== Helper : detect tR =====
function find_tR(sol; smooth_window::Int = 5, rise_frac::Float64 = 0.01)
    ts  = sol.t
    H   = sol[1, :] .+ sol[2, :]

    n   = length(H)
    Hs  = similar(H)
    hw  = smooth_window ÷ 2
    for i in 1:n
        lo = max(1, i - hw); hi = min(n, i + hw)
        Hs[i] = mean(H[lo:hi])
    end

    idx_min = argmin(Hs)
    H_min   = Hs[idx_min]
    H_init  = Hs[1]

    H_min > 0.95 * H_init && return NaN

    threshold = H_min * (1.0 + rise_frac)

    for i in (idx_min + 1):n
        if Hs[i] >= threshold
            t1, t2 = ts[i-1], ts[i]
            H1, H2 = Hs[i-1], Hs[i]
            return t1 + (threshold - H1) / (H2 - H1) * (t2 - t1)
        end
    end
    return NaN
end

function find_slope_before_tR(sol, tR; fit_duration::Float64 = 3.0)
    isnan(tR) && return NaN

    ts = sol.t
    H  = sol[1, :] .+ sol[2, :]

    # Fenêtre : tR - fit_duration → tR (ou fin de simulation)
    mask = (ts .<= tR) .& (ts .>= tR - fit_duration) .& (H .> 0)
    sum(mask) < 3 && return NaN

    t_fit = ts[mask]
    logH  = log.(H[mask])

    # Régression linéaire : logH = a + slope * t
    t_mean   = mean(t_fit)
    logH_mean = mean(logH)
    slope = sum((t_fit .- t_mean) .* (logH .- logH_mean)) /
            sum((t_fit .- t_mean) .^ 2)

    return slope   # unité : jour⁻¹  (taux de croissance net)
end

function find_slope_after_tR(sol, tR; fit_duration::Float64 = 15.0)
    isnan(tR) && return NaN

    ts = sol.t
    H  = sol[1, :] .+ sol[2, :]

    # Fenêtre : tR → tR + fit_duration (ou fin de simulation)
    mask = (ts .>= tR) .& (ts .<= tR + fit_duration) .& (H .> 0)
    sum(mask) < 3 && return NaN

    t_fit = ts[mask]
    logH  = log.(H[mask])

    # Régression linéaire : logH = a + slope * t
    t_mean   = mean(t_fit)
    logH_mean = mean(logH)
    slope = sum((t_fit .- t_mean) .* (logH .- logH_mean)) /
            sum((t_fit .- t_mean) .^ 2)

    return slope   # unité : jour⁻¹  (taux de croissance net)
end

## ===== Test find_tR =====
test_cases = [
    (mop=1e-1, H0=1e6),(mop=1e0, H0=1e6),(mop=5, H0=1e6),(mop=1e1, H0=1e6),(mop=1e2, H0=1e6),(mop=1e3, H0=1e6),
]

n_tc  = length(test_cases)
n_col = 3
n_row = ceil(Int, n_tc / n_col)

pl_test = plot(
    layout        = (n_row, n_col),
    size          = (500 * n_col, 400 * n_row),
    left_margin   = 12mm, right_margin = 5mm,
    top_margin    = 8mm,  bottom_margin = 10mm,
    guidefontsize = 13, tickfontsize = 11, titlefontsize = 12,
    legendfontsize= 10,
)

println("\n── find_tR diagnostic ──")
for (k, tc) in enumerate(test_cases)
    H0  = tc.H0
    mop = tc.mop
    S0  = prop_S_0 * H0
    V0  = mop * H0
    u0  = [S0, 0.0, V0]

    prob = ODEProblem(ModelDef.ODE_MODEL!, u0, t_span, [a, b])
    sol  = solve(prob, Rodas5();
        reltol=1e-6, abstol=1e-6,
        saveat=t_save, isoutofdomain=isoutofdomain, maxiters=1_000_000)

    H     = sol[1, :] .+ sol[2, :]
    tR    = find_tR(sol)
    slope_after = find_slope_after_tR(sol, tR)
    slope_before = find_slope_before_tR(sol, tR)

    H_at_tR = if !isnan(tR)
        idx = clamp(searchsortedfirst(sol.t, tR), 1, length(sol.t))
        H[idx]
    else
        NaN
    end

    @printf "  MOP=%.1e  H0=%.1e  →  tR = %s  |  slope_before = %s  |  slope_after = %s\n" mop H0 (
        isnan(tR)    ? "NaN (no recovery)"      : @sprintf("%.2f d", tR)) (
        isnan(slope_before) ? "NaN"                     : @sprintf("%.4f d⁻¹", slope_before)) (
        isnan(slope_after) ? "NaN"                     : @sprintf("%.4f d⁻¹", slope_after))

    ttl = @sprintf("MOP=%.0e  H₀=%.0e", mop, H0)

    # ── Courbe H(t) ──
    plot!(pl_test[k], sol.t, H,
        label="H(t)", color=:steelblue, lw=2,
        yscale=:log10,
        xlabel="Time [days]", ylabel="H [cells/mL]",
        title=ttl, legend=:bottomright)

    if !isnan(tR)
        vline!(pl_test[k], [tR],
            label=@sprintf("tR=%.1f d", tR),
            color=:red, lw=2, ls=:dash)
        scatter!(pl_test[k], [tR], [H_at_tR],
            label="", color=:red, markersize=7, markerstrokewidth=0)

        # ── Overlay régression linéaire (log-scale) ──
        if !isnan(slope_after)
            fit_duration = 20.0
            t_fit_end    = min(tR + fit_duration, sol.t[end])
            t_reg        = range(tR, t_fit_end, length=50)

            # ancrage : H à tR
            idx_tR  = clamp(searchsortedfirst(sol.t, tR), 1, length(sol.t))
            logH_tR = log(H[idx_tR])
            H_reg   = exp.(logH_tR .+ slope_after .* (t_reg .- tR))

            plot!(pl_test[k], collect(t_reg), H_reg,
                label=@sprintf("slope=%.3f d⁻¹", slope_after),
                color=:orange, lw=2, ls=:dash)
        end
        if !isnan(slope_before)
            fit_duration = 3.0
            t_fit_end    = max(tR - fit_duration, sol.t[1])
            t_reg        = range(t_fit_end, tR, length=50)

            # ancrage : H à tR
            idx_tR  = clamp(searchsortedfirst(sol.t, tR), 1, length(sol.t))
            logH_tR = log(H[idx_tR])
            H_reg   = exp.(logH_tR .+ slope_before .* (t_reg .- tR))

            plot!(pl_test[k], collect(t_reg), H_reg,
                label=@sprintf("slope=%.3f d⁻¹", slope_before),
                color=:orange, lw=2, ls=:dash)
        end
    else
        annotate!(pl_test[k],
            sol.t[end] * 0.5, maximum(filter(isfinite, H)) * 0.5,
            text("no recovery", :red, 10))
    end
end

fig_test = joinpath(output_path, "test_tR_slopes.png")
savefig(pl_test, fig_test)
println("  → Diagnostic figure saved: $fig_test\n")

## ===== 2-D scan =====
res_mop = Float64[]
res_H0  = Float64[]
res_tR  = Float64[]
res_slope_after = Float64[]
res_slope_before = Float64[]

n_total = length(H0_values) * length(mop_values)
global n_done  = 0

println("\nScanning $(length(H0_values)) H0 × $(length(mop_values)) MOP = $n_total simulations …\n")

for H0 in H0_values
    S0 = prop_S_0 * H0
    for mop in mop_values
        V0 = mop * H0
        u0 = [S0, 0.0, V0]

        prob = ODEProblem(ModelDef.ODE_MODEL!, u0, t_span, [a, b])
        sol  = solve(prob, Rodas5();
            reltol=1e-6, abstol=1e-6,
            saveat=t_save,
            isoutofdomain=isoutofdomain,
            maxiters=1_000_000)

        tR = find_tR(sol)
        slope_after = find_slope_after_tR(sol, tR)
        slope_before = find_slope_before_tR(sol, tR)
        push!(res_mop, mop)
        push!(res_H0,  H0)
        push!(res_tR,  tR)
        push!(res_slope_after, slope_after)
        push!(res_slope_before, slope_before)

        global n_done += 1
        if n_done % 20 == 0
            @printf "  [%3d/%3d]  H0=%.1e  MOP=%.2e  tR=%s\n" n_done n_total H0 mop (
                isnan(tR) ? "NaN" : @sprintf("%.1f d", tR))
        end
    end
end

println("\nDone. Recovery detected: $(count(!isnan, res_tR)) / $(length(res_tR))")

## ===== Color palette: one color per H0 =====
using Plots.Colors
palette_colors = cgrad(:viridis, length(H0_values), categorical=true)

## ===== Combined figure: 4 subplots (1 row × 4 cols) =====
xlims_      = (1e-1, 1e10)
xtick_vals  = [10.0^i for i in -1:10]
xtick_labels = [L"10^{%$i}" for i in -1:10]

pl_combined = plot(
    layout        = (1, 4),
    size          = (6000, 1500),
    margins       = 30mm,
    grid          = true,
    titlefontsize  = 30,
    guidefontsize  = 25,
    tickfontsize   = 25,
    legendfontsize = 25,
    xtickfontsize  = 25,
)

# ── Subplot 1 : α = a·MOP + b vs MOP ──
mop_range = 10 .^ range(-4, 2, length=300)
alpha_vals = a .* mop_range .+ b

plot!(pl_combined[1],
    mop_range, alpha_vals,
    xscale  = :log10,
    yscale  = :log10,
    xlabel  = "Initial MOP (V₀ / H₀)",
    ylabel  = "Initial alpha = a.MOP + b",
    title   = "Initial alpha vs MOP",
    label   = "Initial alpha",
    color   = :steelblue,
    lw      = 3,
    legend  = :topleft,
)

# Ligne horizontale b (plancher)
hline!(pl_combined[1], [b],
    label  = "b = $b_str",
    color  = :red,
    lw     = 2,
    ls     = :dash,
)

# Annoter les valeurs de a et b
annotate!(pl_combined[1],
    1e-3, maximum(alpha_vals) * 0.5,
    text("a = $a_str\nb = $b_str", :left, 18, :black),
)

ylabels = [
    "Recovery time  tR  [days]",
    "Growth rate after tR  [day⁻¹]",
    "Growth rate before tR  [day⁻¹]",
]
titles = [
    "Recovery time  tR vs Initial MOP",
    "Growth rate after tR vs Initial MOP",
    "Growth rate before tR vs Initial MOP",
]
legends = [:bottomleft, :bottomleft, :bottomleft]
data_cols = [res_tR, res_slope_after, res_slope_before]

for (sp, (ylab, leg, ydata, title)) in enumerate(zip(ylabels, legends, data_cols, titles))
    plot!(pl_combined[sp + 1],          # décalage de 1 pour laisser la place au subplot α
        xscale  = :log10,
        xlabel  = "Initial MOP (V₀ / H₀)",
        ylabel  = ylab,
        title   = title,
        legend  = leg,
        xlims   = xlims_,
        xticks  = (xtick_vals, xtick_labels),
    )

    for (k, H0) in enumerate(H0_values)
        mask = (res_H0 .== H0) .& (.!isnan.(ydata))
        col  = palette_colors[k]
        e    = floor(Int, log10(H0))
        m    = H0 / 10.0^e
        lbl  = m ≈ 1.0 ? "H₀ = 1e$e" : "H₀ = $(Int(m))e$e"

        scatter!(pl_combined[sp + 1], res_mop[mask], ydata[mask],
            label            = lbl,
            color            = col,
            markersize       = 7,
            markerstrokewidth = 0.5,
            markerstrokecolor = :white,
            alpha            = 0.85,
        )
    end
end

fig_combined = joinpath(output_path, "$(a_str)_$(b_str).png")
savefig(pl_combined, fig_combined)
println("Figure saved → $fig_combined")