## ===== Model : SRVi_phi =====
# dS/dt  = r·S·(1 - H/K) - φ(t)·S·Vi - α·S
# dR/dt  = r·R·(1 - H/K) + α·S
# dVi/dt = β·φ(t)·S·Vi - δ·Vi
#
# State : Y = [S, R, Vi]
# H = S + R
#
# φ(t) : cubic spline in log-space, knots fixed at cycle boundaries
#         log(φ(t_k)) = θ_k  (fitted)
#
# φ_equiv(t) = φ(t) · S(t) / H(t)
# → cohérence avec φ'(t) mesurée expérimentalement

module ModelDef

using DataInterpolations

export MODEL_NAME, FIXED_PARAMS, FITTED_PARAMS, LOG_FIT, LOG_INDICES,
       ODE_MODEL!, PHI_equiv, INITIAL_CONDITION, STATE_LABELS,
       KNOTS, build_phi_func

# ── Identity ────────────────────────────────────────────────────────────────
const MODEL_NAME = "SRVi_phi"

# ── Fixed parameters ────────────────────────────────────────────────────────
const FIXED_PARAMS = (
    r = 0.574619342477644,
    K = 6.675449070379925e7,
    β = 144.0,
    δ = 0.02,
)

# ── Spline knots (cycle boundary times, days) ────────────────────────────────
# One knot per cycle boundary : t_start of cycle 1, t_end of each cycle.
# These must match the time grid used in main.jl (raw_data tH endpoints).
# Edit here if your cycle times differ.
const KNOTS = [0.0000, 6.1927, 12.3854, 18.5781, 24.7708, 27.9931, 31.2153, 34.4375, 37.3264, 40.2153, 43.1042, 46.2917, 49.4792, 52.6667, 55.8542, 59.1181, 62.3819, 65.6458]

const N_KNOTS = length(KNOTS)

# ── Fitted parameters ────────────────────────────────────────────────────────
#
# Layout of the optimiser vector θ_opt (length = 1 + N_KNOTS) :
#   θ_opt[1]          → log(α)          fitted in log-space
#   θ_opt[2:N_KNOTS+1] → log(φ(t_k))   already in log-space, fitted directly
#
# LOG_FIT   = false  : main.jl will NOT apply exp() globally
# LOG_INDICES        : indices where main.jl must apply exp() before passing p
#                      to ODE_MODEL! and PHI_equiv
#                      → only index 1 (α)

const LOG_FIT     = false
const LOG_INDICES = [1]          # only α is stored as log(α) in the optim vector

# FITTED_PARAMS drives bounds display / logs in main.jl
# For θ_k (log-space φ values) bounds are given directly in log-space
const FITTED_PARAMS = vcat(
    [(name = :α, lower = log(1e-30), upper = log(1e-25),
      description = "log mutation rate S→R")],
    [(name = Symbol("logφ_$(k)"), lower = log(1e-14), upper = log(1e-7),
      description = "log φ at knot $(KNOTS[k]) d")
     for k in 1:N_KNOTS]
)

# ── Spline builder ───────────────────────────────────────────────────────────
# p  : full natural-space parameter vector [α, θ_1 … θ_N_KNOTS]
#      (α already decoded by main.jl via LOG_INDICES; θ_k still in log-space)
# Returns a closure φ(t) → Float64
function build_phi_func(p)
    θ = p[2:end]

    spline = CubicSpline(θ, KNOTS)

    t_lo = first(KNOTS)
    t_hi = last(KNOTS)

    return t -> exp(spline(clamp(t, t_lo, t_hi)))
end

# ── ODE ─────────────────────────────────────────────────────────────────────
# p = [α, θ_1 … θ_N_KNOTS]  (α in natural space, θ_k in log-space)
function ODE_MODEL!(dY, Y, p, t)
    α      = p.α
    φt     = p.φ(t)
    S, R, Vi = Y[1], Y[2], Y[3]
    H = S + R
    r, K, β, δ = FIXED_PARAMS.r, FIXED_PARAMS.K, FIXED_PARAMS.β, FIXED_PARAMS.δ
    dY[1] = r*S*(1 - H/K) - φt*S*Vi - α*S
    dY[2] = r*R*(1 - H/K) + α*S
    dY[3] = β*φt*S*Vi - δ*Vi
end

# ── φ_equiv ──────────────────────────────────────────────────────────────────
# Apparent infection rate comparable to φ'(t).
# Y_t : state vector at a single time point
# p   : full parameter vector (same layout as ODE_MODEL!)
# t   : time (needed to evaluate φ(t))
function PHI_equiv(Y_t, p, t)
    φ_func = build_phi_func(p)
    φt = φ_func(t)
    S  = Y_t[1]
    H  = Y_t[1] + Y_t[2]
    return H > 1e-30 ? φt * S / H : φt
end

# ── Initial condition builder ─────────────────────────────────────────────────
function INITIAL_CONDITION(H0, V0, prop_S, prop_Vi)
    return [prop_S * H0, (1 - prop_S) * H0, V0]
end

# ── State labels (for plots / logs) ──────────────────────────────────────────
const STATE_LABELS = [
    (state = :H_total, indices = (1, 2), label = "Model H",  color = :black, lw = 4, ls = :solid),
    (state = :S,       indices = (1,),   label = "Model S",  color = :red,   lw = 2, ls = :dash),
    (state = :R,       indices = (2,),   label = "Model R",  color = :green, lw = 2, ls = :dash),
    (state = :V_total, indices = (3,),   label = "Model V",  color = :black, lw = 4, ls = :solid),
    (state = :Vi,      indices = (3,),   label = "Model Vi", color = :red,   lw = 2, ls = :dash),
]

end # module